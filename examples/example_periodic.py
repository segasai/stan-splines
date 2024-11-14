import cmdstanpy
import matplotlib.pyplot as plt
import numpy as np
import scipy.interpolate

np.random.seed(11)

N = 100
x = np.random.uniform(0, 2 * np.pi, size=N)
y0 = 1 + np.sin(x) * 3
ey = x * 0 + 0.3
y = y0 + np.random.normal(size=N) * ey
xknots = np.array([0, 0.5, 1, 2, 4, 2 * np.pi])
nknots = len(xknots)
#M = stan.build(open('example.stan').read())
M = cmdstanpy.CmdStanModel('example_periodic', 'example_periodic.stan')
data = {'N': N, 'x': x, 'y': y, 'ey': ey, 'xknots': xknots, 'nknots': nknots}
R = M.sample(data=data, seed=434)
res = R.stan_variables()
plt.plot(x, y, '.')
nplot = 30
xgrid = np.linspace(0, 2 * np.pi, 1000)
means = res['yknots'].mean(axis=0)
stds = res['yknots'].std(axis=0)
means = np.concatenate((means, [means[0]]))
stds = np.concatenate((stds, [stds[0]]))
plt.errorbar(xknots, means, stds, fmt='.', color='black')
for i in range(nplot):
    C = scipy.interpolate.CubicSpline(xknots,
                                      np.concatenate((res['yknots'][i],
                                                      [res['yknots'][i][0]])),
                                      bc_type='periodic')
    plt.plot(xgrid, C(xgrid), alpha=0.2, color='red')

plt.savefig('plot_example_periodic.png')
