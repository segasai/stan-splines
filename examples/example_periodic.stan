functions
{
#include ../stan/spline_periodic.stan
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
  array[N] int x_pos_knots = spline_findpos(xknots, x);
}
parameters
{
  // the parameters of our spline model are
  // the values at the knots
  vector[nknots-1] yknots;
}
transformed parameters
{
  vector[nknots] spl_coeffs = spline_getcoeffs(xknots, yknots);
  // these are the spline coefficients corresponding to the current model
}

model
{
  vector[N] ymod;
  ymod = spline_eval(xknots,
		     yknots, spl_coeffs, x, x_pos_knots);
  y ~ normal (ymod, ey);
}
generated quantities
{
  vector[N] ymod;
  ymod = spline_eval(xknots,
		     yknots, spl_coeffs, x, x_pos_knots);
}
