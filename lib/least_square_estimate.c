#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_fit.h>
#include <gsl/gsl_cdf.h>

int
least_square_estimate(const double *x, const double *y, size_t n, double *res) {
    // y = c0 + c1 * x
    // return c1, stdev1, t1, p_value
    double c0, c1, cov00, cov01, cov11, sumsq;
    gsl_fit_linear(x, 1, y, 1, n, &c0, &c1, &cov00, &cov01, &cov11, &sumsq);
    printf("%lf %lf %lf %lf %lf %lf\n", c0, c1, cov00, cov01, cov11, sumsq);
    
    //double stdev0 = sqrt(cov00);
    //double t0 = c0 / stdev0;
    //double pv0 = (t0 < 0)? 2 * (1 - gsl_cdf_tdist_P(-t0, n - 2)): 2 * (1 - gsl_cdf_tdist_P(t0, n - 2));

    double stdev1 = sqrt(cov11);
    double t1 = c1 / stdev1;
    //double pv1 = t1 < 0 ? 2 * (1 - gsl_cdf_tdist_P(-t1, n - 2)): 2 * (1 - gsl_cdf_tdist_P(t1, n - 2));

    int i = 0;
    double dl = n - 2;
    double y_mean = 0;
    for (i = 0; i < n; i++) {
        y_mean += y[i];
    }
    y_mean = y_mean / n;
    
    double y_var = 0;
    for (i = 0; i < n; i++) {
        y_var +=  pow(y[i] - y_mean, 2);
    }
    
    //double ym = 0.2 * (y[0] + y[1] + y[2] + y[3] + y[4]);
    //double sct = pow(y[0] - ym, 2) + pow(y[1]-ym, 2) + pow(y[2] - ym, 2) + pow(y[3] - ym, 2) + pow(y[4] - ym, 2);
    double R2 = 1 - sumsq / y_var;
    double F = R2 * dl / (1 - R2);
    double p_value = 1 - gsl_cdf_fdist_P(F, 1, dl);

    printf("%lf %lf %lf %lf\n", c1, stdev1, t1, p_value);
    res[0] = c1;
    res[1] = stdev1;
    res[2] = t1;
    res[3] = p_value;

    return 0;
}

#define TEST
#ifdef TEST
int
main(void)
{
    double x[5] = {1., 2., 3., 4., 5.};
    double y[5] = {1.1, 1.9, 3.1, 4.2, 4.9};
    size_t n = 5;
    double res[4];
    least_square_estimate(x, y, n, res);

}
#endif
