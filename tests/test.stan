#include spline.stan

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
	int x_pos_knots[N] = findpos(nknots, xknots, N, x);
}
parameters {}
transformed parameters
{
	vector[nknots] spl_coeffs = getcoeffs(nknots, xknots, yknots);
	// these are the spline coefficients corresponding to the current model
}

model
{
  vector[N] ymod;
  ymod = spline_eval(nknots, xknots,
		     yknots, spl_coeffs, N, x, x_pos_knots);
  for (i in 1:N)
    {
      print(ymod[i]);
    }
}