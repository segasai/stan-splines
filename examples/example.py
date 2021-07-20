import cmdstanpy
import matplotlib.pyplot as plt
import numpy as np
import scipy.interpolate

print('Please make sure to copy the spline.stan in this directory!!!')
np.random.seed(11)

N = 100
x = np.random.uniform(0, 1, size=N)
y0 = 1 + x + 4 * x**2 + np.sin(x * 3) * 3
ey = x * 0 + 0.3
y = y0 + np.random.normal(size=N) * ey
xknots = np.array([0, 0.2, 0.4, 0.7, .9, 1])
nknots = len(xknots)
#M = stan.build(open('example.stan').read())
M = cmdstanpy.CmdStanModel('example', 'example.stan')
data = {'N': N, 'x': x, 'y': y, 'ey': ey, 'xknots': xknots, 'nknots': nknots}
R = M.sample(data=data, seed=434)
res = R.stan_variables()
plt.plot(x, y, '.')
nplot = 30
xgrid = np.linspace(0, 1, 1000)
plt.errorbar(xknots,
             res['yknots'].mean(axis=0),
             res['yknots'].std(axis=0),
             fmt='.',
             color='black')
for i in range(nplot):
    C = scipy.interpolate.CubicSpline(xknots,
                                      res['yknots'][i],
                                      bc_type='natural')
    plt.plot(xgrid, C(xgrid), alpha=0.2, color='red')

plt.savefig('result.png')
