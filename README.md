# stan-splines
This is the scode written to implement natural spline models in STAN.

Author: Sergey Koposov

Email: skoposov AT ed DOT ac DOT uk

Please cite this https://ui.adsabs.harvard.edu/abs/2019MNRAS.485.4726K/abstract
paper if you use it.


The splines are definite by the location of the knots and values there.

The usage should be pretty simple. Here is an example of the model that
fits the set of x,y by a spline.

```stan
#include spline.stan

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
	int x_pos_knots[N] = findpos(nknots, xknots, N, x);
}
parameters
{
	// the parameters of our spline model are
	// the values at the knots
	vector[nknots] yknots;
}
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
  y ~ normal (ymod, ey);
}
```

See the jupyter notebook here 
https://github.com/segasai/stan-splines/blob/main/examples/Example.ipynb
