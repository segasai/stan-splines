import cmdstanpy
import matplotlib.pyplot as plt
import numpy as np
import scipy.interpolate

np.random.seed(11)

nknots = 10
N = 1000
xknots = np.concatenate(
    ([0], np.sort(np.random.uniform(0, 2 * np.pi,
                                    size=nknots - 2)), [2 * np.pi]))
yknots = np.random.normal(size=nknots - 1)
x = np.random.uniform(0, 2 * np.pi, size=N)
M = cmdstanpy.CmdStanModel(stan_file='test_periodic.stan')
data = {'xknots': xknots, 'nknots': nknots, 'x': x, 'N': N, 'yknots': yknots}
R = M.optimize(data=data, seed=434, iter=1)
fout = R.runset.stdout_files[0]
with open(fout, 'r') as fp:
    header = True
    res = []
    while True:
        curl = fp.readline().rstrip()
        if header:
            try:
                val = float(curl)
            except ValueError:
                continue
            header = False
            res.append(val)
        else:
            try:
                val = float(curl)
                res.append(val)
            except ValueError:
                break
res = np.array(res)
C = scipy.interpolate.CubicSpline(xknots,
                                  np.concatenate([yknots, [yknots[0]]]),
                                  bc_type='periodic')
ypred = C(x)
maxdiff = np.abs(ypred - res).max()
print('MAXIMUM difference is ', maxdiff)
assert (maxdiff < 1e-5)
