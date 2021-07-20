#include spline_precompute.stan

data{
  int nknots;
  int N;
  vector[nknots] xknots;
  vector[nknots] yknots;
  vector[N] x;
}
transformed data
{
  // determine which knots the point belong to
  int x_pos_knots[N] = spline_findpos(nknots, xknots, N, x);
  // precomuted the polynomials
  vector[N] spline_mat[4] = spline_getmat(x, N, xknots, nknots, x_pos_knots);
  
}
parameters {}
transformed parameters
{
  vector[nknots] spl_coeffs = spline_getcoeffs(nknots, xknots, yknots);
  // these are the spline coefficients corresponding to the current model
}

model
{
  vector[N] ymod;
  ymod = spline_eval(nknots, xknots,
		     yknots, spl_coeffs, N, spline_mat, x_pos_knots);
  for (i in 1:N)
    {
      print(ymod[i]);
    }
}
