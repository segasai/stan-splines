functions
{
#include ../stan/spline.stan
}
data{
  int N;
  int nknots;
  vector[N] x;
  vector[N] y;
  vector[N] ey;
  vector[nknots] xknots;
}
transformed data
{
  // determine which knots the point belong to
  int x_pos_knots[N] = spline_findpos(nknots, xknots, N, x);
}
parameters
{
  // the parameters of our spline model are
  // the values at the knots
  vector[nknots] yknots;
}
transformed parameters
{
  vector[nknots] spl_coeffs = spline_getcoeffs(nknots, xknots, yknots);
  // these are the spline coefficients corresponding to the current model
}

model
{
  vector[N] ymod;
  ymod = spline_eval(nknots, xknots,
		     yknots, spl_coeffs, N, x, x_pos_knots);
  y ~ normal (ymod, ey);
}
