/*Nonlinear Least Squares Curve Fitting Program*/

/*Marquardt algorithm from P.R. Bevington,"Data Reduction and Error
Analysis for the Physical Sciences," McGraw-Hill, 1969; Adapted by
Wayne Weimer & David Harris. Jackknife error algorithm of M.S. Caceci,
Anal. Chem. 1989, 61, 2324. Translated to ANSI C by Douglas Harris &
Tim Seufert 7-94 */
#include "stdafx.h"
#include <stdio.h>
#include <Stdlib.h>
#include <math.h>
#include <ctype.h>

#define _CRT_SECURE_NO_WARNINGS
#define maxnpts 50 /* Maximum data pairs -increase if desired */

/*Change nterms to the number of parameters to be fit in your equation*/
/***********************************************************/
#define nterms 2
/* Number of parameters to be fit */
/***********************************************************/
int param, iteration, nloops, n, cycle, nfree;
int npts;                                                   /* Number of data pairs */
double x[maxnpts], y[maxnpts], sigmay[maxnpts];        /*x,y,y uncertainty*/
double weight[maxnpts];                                /*Weighting factor*/
double yfit[maxnpts];                                  /*Calculated values of y */
double a[nterms];                                      /* a[i]=c[i] params */
double sigmaa[nterms];
double b[nterms];
double beta[nterms], c[nterms];                        /*To be fit by program*/
double finala[nterms], lastsigmaa[nterms];
double alpha[nterms][nterms];
double arry[nterms][nterms];
double aug[nterms][nterms * 2];                        /* For matrix inversion */
double deriv[maxnpts][nterms];                         /* Derivatives */
double flambda;                                        /*Proportion of gradient search(=0.001 at start)*/
double chisq;                                          /* Variance of residuals in curve fit */
double chisql, fchisq, sy;
char errorchoice;
const int BUFF_SIZE = 100;
char filename[20];
char answer[BUFF_SIZE];
FILE *fp;
void readdata();
void unweightedinput();
void weightedinput();
void chisquare();
void calcderivative();
void matrixinvert();
void curvefit(int npoints);
void display();
void uncertainties();
void jackknifedata(char *filename, int k);
void print_matrix(double matirx[][nterms], int size_y);
void print_array(double _arry[], int size);

#if defined _WIN32
errno_t err;
#endif

int main() {
    int i;
    printf("Least Squares Curve Fitting. You must modify the constant\n");
    printf("'nterms' and the fuction 'Fune' for new problems.\n");
    readdata();
    printf("\nEnter initial guesses for parameters:\n");
    printf("\t(Note: Parameters cannot be exactly zero.)\n");
    for (i = 0; i < nterms; i++) {
        while(a[i] == 0.0) {
            printf("Parameter #%d =   ", i + 1);
            fgets(answer, BUFF_SIZE, stdin);
            a[i] = atof(answer);
        }
    }
    printf("Initial A array:\n");
    print_array(a, nterms);
    flambda = 0.001;
    iteration = 0;
    cycle = 0;
    do {
        curvefit(npts);
        iteration++;
        display();
        iteration = 0;
        cycle = 0;
        printf("\n\tAnother iteration (Y/N)? ");
        fgets(answer, BUFF_SIZE, stdin);
    } while (answer[0] != 'N' && answer[0] != 'n');
    printf("\nDo you want to calculate uncertainty in parameters (Y/N)?");
    fgets(answer, BUFF_SIZE, stdin);
    if (answer[0] == 'Y' || answer[0] == 'y') uncertainties();
    return 0;
}

// Displays the data entered
void print_data() {
    int i;
    for (i = 0; i < npts; i++) {
        printf("%d\tx = %- #12.8f\ty = %- #12.8f\t", i + 1, x[i], y[i]);
        printf("Sigmay = %- #12.8f\n", sigmay[i]);
        weight[i] = 1 / ( sigmay[i] * sigmay[i] );
    }
}

                        /*******************************/
double func(int i) /* The function you are fitting*/
{                       /*******************************/
    int loop;
    double value;
    if (param == 1) {
        for (loop = 0; loop < nterms; loop++) {
            c[loop] = b[loop];
        }
    }
    else {
        for (loop = 0; loop < nterms; loop++) {
            c[loop] = a[loop];
        }
    }

    /********************************************/
    /*      Enter the function to be fit:       */
    /********************************************/    
    value = c[0] * x[i] + c[1]; /*Linear Equation*/    
    printf("\nfunc(i) -- i: %d  a: %f  x: %f  b: %f  =  value: %f\n", i, c[0], x[i], c[1], value);
    // x[i] is the independent variable
    // Values of c[n], c[n-1], c[0] are determined by least squares
    // nterms must be set to equal n+l
    // Example of another function: value= c[2]*x[i]*x[i]+c[l]*x[i]+c[O]
    return ( value );
}

void readdata() {
    int n = 0;

    // Prompt for data entry type
    do {
        printf("\nDo you want to enter x,y values or read them from a file?\n");
        printf("\tType E for enter and F for File: ");
        fgets(answer, BUFF_SIZE, stdin);
        answer[0] = toupper(answer[0]);
    } while (answer[0] != 'E' && answer[0] != 'F');

    // Read from file
    if (answer[0] == 'F') {
        do {
            printf("\nPlease enter the name of the data file: ");
            fgets(filename, BUFF_SIZE, stdin);
            printf("\n");
#if defined _WIN32
            err = fopen_s(&fp, filename, "rb");
            if (err != 0) {
                printf("Fatal error: could not open file %s\n", filename);
                exit(1);
            }
#else
            fp = fopen(filename, "rb");
            if (fp == NULL) {
                printf("Fatal error: could not open file %s\n", filename);
                exit(1);
            }
#endif

            for (n = 0; !feof(fp); n++) {
                fread(&x[n], sizeof( double ), 1, fp);
                fread(&y[n], sizeof( double ), 1, fp);
                fread(&sigmay[n], sizeof( double ), 1, fp);
                if (errorchoice == '1') {
                    sigmay[n] = 1.0;
                }
            }
            fclose(fp);
            npts = n - 1;
            print_data();
            printf("\nIs this data correct (Y/N)?");
            fgets(answer, BUFF_SIZE, stdin);
        } while (answer[0] != 'Y' && answer[0] != 'y');
    }
    // Enter data manually
    else {
        do {
            printf("\nChoices for error analysis : \n");
            printf("\tl. Let the program weight all points equally\n");
            printf("\t2. Enter estimated uncertainty for each point\n\n");
            printf("Choose option 1 or 2 now: ");
            fgets(answer, BUFF_SIZE, stdin);
        } while (answer[0] != '1' && answer[0] != '2');

        errorchoice = answer[0];

        do {
            if (errorchoice == '1') {
                printf("UW input\n");
                unweightedinput();
            }
            else if (errorchoice == '2') {
                printf("Weighted input\n");
                weightedinput();
            }
            print_data();
            printf("Is this data correct(Y/N)?");
            fgets(answer, BUFF_SIZE, stdin);
        } while (answer[0] != 'y' && answer[0] != 'Y');
        printf("Enter name of file to save the data in: ");
        fgets(filename, BUFF_SIZE, stdin);
#if defined _WIN32
        err = fopen_s(&fp, filename, "wb");
        if (err != 0) {
            printf("Fatal error: could not open file %s\n", filename);
            exit(1);
        }
#else
        fp = fopen(filename, "wb");
        if (fp == NULL) {
            printf("Fatal error: could not open file %s\n", filename);
            exit(1);
        }
#endif

        for (n = 0; n < npts; n++) {
            fwrite(&x[n], sizeof( double ), 1, fp);
            fwrite(&y[n], sizeof( double ), 1, fp);
            fwrite(&sigmay[n], sizeof( double ), 1, fp);
        }

        fclose(fp);
        printf("Data saved in file %s\n", filename);
    }
}

/* Enter equal weight data */
void unweightedinput() {
    int i, n;
    printf("List the data in the order: x y, with one set on each\n");
    printf("line and a space (not a comma) between the numbers.\n");
    printf("Type END to end input\n");
    for (n = 0;; n++) {
        fgets(answer, BUFF_SIZE, stdin);
        if (answer[0] == 'E' || answer[0] == 'e') {
            break;
        }
        // Convert first part of string input
        x[n] = atof(answer);
        i = 0;
        while (answer[i] != ' ' && answer[i] != '\0') {
            i++;
        }
        // Convert second half of string input
        y[n] = atof(answer + i);
        sigmay[n] = 1;
    }
    npts = n;
}

// Enter unequal weighted data
void weightedinput() {
    int i, n;
    printf("List the data in the order: x   y sigmay, with one set on\n");
    printf("each line and a space (not a comma) between the numbers.\n");
    printf("Type END to end input\n"); for (n = 0;; n++) {
        fgets(answer, BUFF_SIZE, stdin);
        if (answer[0] == 'E' || answer[0] == 'e') {
            break;
        }
        x[n] = atof(answer);
        i = 0;
        while (answer[i] != ' ' && answer[i] != '\0') {
            i++;
        }
        y[n] = atof(answer + i);
        i++;
        while (answer[i] != ' ' && answer[i] != '\0') {
            i++;
        }
        sigmay[n] = atof(answer + i);
    }
    npts = n;
}

// Sum of square of differences between measured and calculated y values
void chisquare() {
    int i;
    fchisq = 0;
    for (i = 0; i < npts; i++){        
        fchisq += weight[i] * ( y[i] - yfit[i] ) * ( y[i] - yfit[i] );
        printf("y[i]: %f -- yfit[i]: %f -- fchisq: %f\n", y[i], yfit[i], fchisq);
    }
    fchisq /= nfree;
    printf("Final chisq: %f\n\n", fchisq);
}

// Numerical derivative
void calcderivative() {
    int i, m;
    double atemp, delta;
    for (m = 0; m < nterms; m++) {
        atemp = a[m];
        delta = fabs(a[m] / 100000);
        a[m] = atemp + delta;
        for (i = 0; i < npts; i++) {
            deriv[i][m] = ( func(i) - yfit[i] ) / delta;
            a[m] = atemp;
        }
        printf("\nNumerical derivative matrix:\n");
        print_matrix(deriv, npts);
    }
}

// Inverts the matrix array[][]
// Pivoting reduces rounding error
void matrixinvert() {
    int i, j, k, ik[nterms], jk[nterms];
    double rsave, amax;

    for (k = 0; k < nterms; k++) {

        amax = 0.0;

        for (i = k; i < nterms; i++) {
            for (j = k; j < nterms; j++) {
                if (abs(amax) <= abs(arry[i][j])) {
                    amax = arry[i][j];
                    ik[k] = i;
                    jk[k] = j;
                }
            }
        }

        i = ik[k];

        if (i > k) {
            for (j = 0; j < nterms; j++) {
                rsave = arry[k][j];
                arry[k][j] = arry[i][j];
                arry[i][j] = -1 * rsave;
            }
        }

        j = jk[k];

        if (j>k) {
            for (i = 0; i < nterms; i++) {
                rsave = arry[i][k];
                arry[i][k] = arry[i][j];
                arry[i][j] = -1 * rsave;
            }
        }
        for (i = 0; i < nterms; i++) {
            if (i != k) {
                arry[i][k] = -1 * arry[i][k] / amax;
            }
        }
        for (i = 0; i < nterms; i++) {
            for (j = 0; j < nterms; j++) {
                if (j != k && i != k) {
                    arry[i][j] = arry[i][j] + arry[i][k] * arry[k][j];
                }
            }
        }
        for (j = 0; j < nterms; j++) {
            if (j != k) {
                arry[k][j] = arry[k][j] / amax;
            }
        }
        arry[k][k] = 1 / amax;
    }
    for (k = nterms - 1; k > -1; k--) {
        j = ik[k];
        if (j > k) {
            for (i = 0; i < nterms; i++) {
                rsave = arry[i][k];
                arry[i][k] = -1 * arry[i][j];
                arry[i][j] = rsave;
            }
        }
        i = jk[k];
        if (i > k) {
            for (j = 0; j < nterms; j++) {
                rsave = arry[k][j];
                arry[k][j] = -1 * arry[i][j];
                arry[i][j] = rsave;
            }
        }
    }
}

// Curve fitting algorithm
void curvefit(int npoints) {
    int i, j, k;
    nfree = npoints - nterms;
    
    // Clear b and beta arrays
    for (j = 0; j < nterms; j++) {
        b[j] = beta[j] = 0;
        for (k = 0; k <= j; k++) {
            alpha[j][k] = 0;
        }
    }
    param = 0;
    
    // Find y values for current parameter values
    for (i = 0; i < npoints; i++) {
        yfit[i] = func(i);
    }
    printf("\nyfit array:\n");
    print_array(yfit, npts);

    // Find the chi squared value
    chisquare();
    chisql = fchisq;

    // Find the derivative
    calcderivative();

    // For each data point set...
    for (i = 0; i < npoints; i++) {
        // ... for each parmeter term...
        for (j = 0; j < nterms; j++) {
            // beta = weight * (data point y - estimated y) * derivative
            beta[j] += weight[i] * ( y[i] - yfit[i] ) * deriv[i][j];
            for (k = 0; k <= j; k++) {
                alpha[j][k] += ( weight[i] * deriv[i][j] * deriv[i][k] );
            }
        }
    }
    printf("\nPopulated alpha array:\n");
    print_matrix(alpha, nterms);
    for (j = 0; j < nterms; j++) {
        for (k = 0; k <= j; k++) {
            alpha[k][j] = alpha[j][k];
        }
    }
    nloops = 0;
    do {
        param = 1;
        for (j = 0; j < nterms; j++) {
            for (k = 0; k < nterms; k++) {
                arry[j][k] = alpha[j][k] / sqrt(alpha[j][j] * alpha[k][k]);
            }
            arry[j][j] = flambda + 1;
        }
        matrixinvert();
        for (j = 0; j < nterms; j++) {
            b[j] += beta[k] * arry[j][k] / sqrt(alpha[j][j] * alpha[k][k]);
        }
        for (i = 0; i < npoints; i++) {
            yfit[i] = func(i);
        }
        chisquare();
        if (( chisql - fchisq ) < 0) {
            flambda *= 10;
        }
        nloops++;
    } while (fchisq > chisql);
    for (j = 0; j < nterms; j++) {
        a[j] = b[j];
        sigmaa[j] = sqrt(arry[j][j] / alpha[j][j]);
    }
    flambda /= 10;
}

// Prints result of curve fit
void display() {
    int i;
    printf("\nIteration #%d\n", iteration);
    for (i = 0; i < nterms; i++) {
        printf("A[%3dl = %-#12.8f\n", i, a[i]);
        finala[i] = a[i];
    }
    printf("Sum of squares of residuals = %- #12.8f", fchisq * nfree);
    sy = sqrt(fchisq);
}

// Calculates uncertainties by removing one data point and recalculating parameters
void uncertainties() {
    int i, k;
    double ajack[nterms][maxnpts];
    double avajack[nterms];

    do {
        cycle++;
        printf("Calculating uncertainties...");
        for (i = 0; i < npts; i++) {
            jackknifedata(filename, i++);
            for (k = 0; k <= iteration; k++) {
                curvefit(npts - 1);
            }
            printf("Now playing with the data point %d\n", i + 1);
            for (k = 0; k < nterms; k++) {
                ajack[k][i] = a[k];
            }
        }
        printf("\n\n");
        for (k = 0; k < nterms; k++) {
            avajack[k] = 0;
        }
        for (k = 0; k < nterms; k++) {
            for (i = 0; i < npts; i++) {
                avajack[k] += ajack[k][i];
            }
            avajack[k] = avajack[k] / npts;
        }
        for (k = 0; k < nterms; k++) {
            sigmaa[k] = 0;
        }
        for (k = 0; k < nterms; k++) {
            for (i = 0; i < npts; i++) {
                sigmaa[k] += ( ajack[k][i] - avajack[k] ) * ( ajack[k][i] - avajack[k] );
            }
            sigmaa[k] = sqrt(( npts - 1 ) * sigmaa[k] / npts);
            printf("Parameter[%d] = %- 12.8f +/- %- 12.8f\n", k, finala[k], sigmaa[k]);
            if (cycle > 1) {
                printf("\t(Previous uncertainty = %- #12.8f)\n\n", lastsigmaa[k]);
                lastsigmaa[k] = sigmaa[k];
            }
        }
        printf("Standard deviation of y = %-#12.8f\n", sy);
        printf("Above result is based %d iterations\n", iteration);
        iteration += 5;
        printf("Iterations will now be increased to %d" " to see if the estimates of \n", iteration);
        printf("uncertainty change. When consecutive cycles give\n");
        printf("similar results, it is time to stop.\n");
        printf("\tDo you want to try another cycle now (Y/N)? ");
        fgets(answer, BUFF_SIZE, stdin);
    } while (answer[0] == 'y' || answer[0] == 'Y');
}

// Removes one data point
void jackknifedata(char *filename, int k) {
    int n = 0;    
#if defined _WIN32
    err = fopen_s(&fp, filename, "rb");
    if (err != 0) {
        printf("Fatal error: could not open file %s\n", filename);
        exit(1);
    }
#else
    fp = fopen(filename, "rb");
    if (fp == NULL) {
        printf("Fatal error: could not open file %s\n", filename);
        exit(1);
    }
#endif

    while (!feof(fp)) {
        fread(&x[n], sizeof( double ), 1, fp);
        fread(&y[n], sizeof( double ), 1, fp);
        fread(&sigmay[n], sizeof( double ), 1, fp);
        if (errorchoice == 'l') {
            sigmay[n] = 1.0;
        }
        weight[n] = 1 / ( sigmay[n] * sigmay[n] );
        n++;
        npts = n - 1;
        fclose(fp);
        for (n = 0; n < ( npts - 1 ); n++) {
            if (n >= k) {
                x[n] = x[n + 1];
                y[n] = y[n + 1];
                weight[n] = weight[n + 1];
            }
        }
    }
}

void print_matrix(double matrix[][nterms], int size_x) {
    for (int i = 0; i < nterms; i++) {
        for (int j = 0; j < size_x; j++) {
            printf("%f, ", matrix[i][j]);
        }
        printf("\n");
    }
}

void print_array(double this_array[], int size) {
    for (int i = 0; i < size; i++) {
        printf("%f, ", this_array[i]);
    }
    printf("\n");
}
