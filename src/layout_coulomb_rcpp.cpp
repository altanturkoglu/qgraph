#include <R.h>
#include <Rcpp.h>
#include <Rmath.h>

using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
NumericMatrix layout_with_coulomb_Cpp(
    int pniter,
    int pvcount,
    int pecount,
    NumericVector maxdelta,
    double parea,
    double pcoolexp,
    double prepulserad,
    double pcoulombrad,
    IntegerVector Ef, /* Edges from */
    IntegerVector Et, /* Edges to */
    NumericVector W,
    NumericVector sizes,
    NumericVector xInit,
    NumericVector yInit,
    LogicalVector Cx,
    LogicalVector Cy,
    int digits
) {
    /*
     Calculate a two-dimensional Fruchterman-Reingold layout for (symmetrized)
     edgelist matrix d.  Positions (stored in (x,y)) should be initialized
     prior to calling this routine.
     */

    int n = pvcount;
    int m = pecount;
    double frk;
    double ded;
    double xd;
    double yd;
    double rf;
    double af;
    int i;
    int j;
    int k;
    int l;
    int niter = pniter;
    // double maxdelta;
    double area = parea;
    double xmax;
    double ymax;
    double coolexp = pcoolexp;
    double repulserad = prepulserad;
    double coulombrad = pcoulombrad;

    /* Allocate memory for transient structures */
    // dx=(double *)R_alloc(n,sizeof(double));
    // dy=(double *)R_alloc(n,sizeof(double));
    // t=(double *)R_alloc(n,sizeof(double));
    // Rcpp way:
    NumericVector dx(n);
    NumericVector dy(n);
    NumericVector t(n);

    // Copy xIint and yInit:
    NumericVector x(n);
    NumericVector y(n);
    for (i = 0; i < n; i++) {
        x[i] = xInit[i];
        y[i] = yInit[i];
    }

    frk = sqrt(area / (double)n);

    /* Run the annealing loop */
    for (i = niter; i >= 0; i--) {

        /* Check for user interrupt */
        R_CheckUserInterrupt();

        /* Clear the deltas */
        for (j = 0; j < n; j++) {
            dx[j] = 0.0;
            dy[j] = 0.0;
        }

        /* Increment deltas for each undirected pair */
        for (j = 0; j < n; j++) {

            /* Set the temperature (maximum move/iteration) */
            t[j] = maxdelta[j] * pow((double)i / (double)niter, coolexp);  // cooling function is (i/niter)^coolexp

            if (j < n) {
                for (k = j + 1; k < n; k++) {

                    /* Obtain difference vector */
                    xd = x[j] - x[k];
                    yd = y[j] - y[k];

                    /* Get dyadic euclidean distance */
                    ded = sqrt(xd * xd + yd * yd);

                    /* Rescale differences to length 1 */
                    xd /= ded;
                    yd /= ded;

                    /* Calculate base repulsive "spring force" */
                    rf = (frk * frk) * (1.0 / ded - ded * ded / repulserad);

                    /* Scale base repulsive force by a "Coulomb factor" */
                    rf *= (sizes[j] * sizes[k]) / (coulombrad);

                    /* Add to the position change vector */
                    dx[j] += xd * rf;
                    dx[k] -= xd * rf;
                    dy[j] += yd * rf;
                    dy[k] -= yd * rf;
                }
            }
        }

        for (j = 0; j < m; j++) {

            k = Ef[j];
            l = Et[j];

            xd = x[k] - x[l];
            yd = y[k] - y[l];

            /* Get dyadic euclidean distance */
            ded = sqrt(xd * xd + yd * yd); 

            /* Rescale differences to length 1 */
            if (abs(ded) > 0.000001) {
                xd /= ded; 
                yd /= ded;
            }

            /* Calculate the attractive "force" */
            af = ded * ded / frk * W[j];
            
            /* Add to the position change vector */
            dx[k] -= xd * af; 
            dx[l] += xd * af;
            dy[k] -= yd * af;
            dy[l] += yd * af;
        }

        /* Dampen motion, if needed, and move the points */
        for (j = 0; j < n; j++) {
            ded = sqrt(dx[j] * dx[j] + dy[j] * dy[j]);

            /* Dampen to t if necessary */
            if (ded > t[j]) {
                ded = t[j] / ded;
                dx[j] *= ded;
                dy[j] *= ded;
            }

            /* Update positions (correcting for floating point errors) */
            if (!Cx[j]) {
                x[j] += Rf_fround(dx[j], digits);
            }
            if (!Cy[j]) {
                y[j] += Rf_fround(dy[j], digits);
            }
        }
    }

    NumericMatrix Layout(n, 2);

    // Fill layout:
    for (i = 0; i < n; i++) {
        Layout(i, 0) = x[i];
        Layout(i, 1) = y[i];
    }

    return Layout;
}