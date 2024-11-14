functions
{
#include ../stan/spline_periodic.stan
}

data{
  int nknots;
  int N;
  vector[nknots] xknots;
  vector[nknots-1] yknots;
  vector[N] x;
}
transformed data
{
  // determine which knots the point belong to
  array[N] int x_pos_knots = spline_findpos(xknots, x);
}
parameters {}
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
  for (i in 1:N)
    {
      print(ymod[i]);
    }
}
