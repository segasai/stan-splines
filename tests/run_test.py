import cmdstanpy
import matplotlib.pyplot as plt
import numpy as np
import scipy.interpolate

print('Please make sure to copy the spline.stan in this directory!!!')
np.random.seed(11)

nknots = 10
N = 1000
xknots = np.concatenate(
    ([0], np.sort(np.random.uniform(size=(nknots - 2))), [1]))
yknots = np.random.normal(size=nknots)
x = np.random.uniform(size=N)
M = cmdstanpy.CmdStanModel('test', 'test.stan')
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
C = scipy.interpolate.CubicSpline(xknots, yknots, bc_type='natural')
ypred = C(x)
maxdiff = np.abs(ypred - res).max()
print('MAXIMUM difference is ', maxdiff)
assert (maxdiff < 1e-5)
