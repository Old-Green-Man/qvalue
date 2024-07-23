import sys
import numpy as np
import matplotlib.pyplot as plt
# import scipy
# from scipy.interpolate import UnivariateSpline

def qvalue(pv, pi0=None, m=None, verbose=False, plot=False):
  """
  NOTE: Other python versions of this available on the web do not properly
  estimate pi0 because they use splrep with s=0 (no smoothing). The pi0 estimate
  from these functions is simply the last pi0 value considered (which is usually
  a poor estimate) and the spline fitting does nothing. e.g.

  https://pypi.org/project/qvalue/
  https://github.com/nfusi/qvalue
  https://gist.github.com/ryananeff/0232597b04ec1e5947de2ad8b9292d6e
  
  Here, instead of splines, we are calculating the standard deviations of the
  pi0s > lambda. Since the pi0s tend to flatten out as lambda approaches 1 the
  standard deviation tends to be at a minimum at the lambda where the pi0s
  converge to the null proportion. The pi0 is taken as the mean of the pi0s >
  lambda where the standard deviation is at a minimum. See
  https://github.com/Old-Green-Man/qvalue/blob/main/README.md for a graphical
  illustration.

  It is probably a good idea to call this with plot=True on your p-values at
  least once, to see what the distribution looks like. If for some reason it
  seems off you can try providing your own pi0 with the pi0=my_pi0 option. A higher pi0
  will give more conservative q-value estimates.

  Estimates q-values from p-values
  Args
  =====
  pi0: if None, it's estimated based on principles discussed in Storey and Tibshirani, 2003.
  m: number of tests. If not specified m = pv.size
  verbose: print verbose messages? (default False)

  returns 1. a numpy array of q-values with the same shape and order as the input p-values, and 2. the pi0 estimate
  """
  assert(pv.min() >= 0 and pv.max() <= 1), "p-values should be between 0 and 1"

  original_shape = pv.shape
  pv = pv.ravel()  # flattens the array in place, more efficient than flatten()

  if m is None:
    m = float(len(pv))
  else:
    # the user has supplied an m
    m *= 1.0

  pi0 = estimate_pi0(pv, m, pi0, plot=plot)
  # for pwr in range(-10, 11):
  #   Note: s=None is used in multipy
  #   pi0 = multipy_est_pi0(pv, m, s=2**pwr, plot=plot)
  

  p_ordered = np.argsort(pv)
  pv = pv[p_ordered]
  qv = pi0 * m / len(pv) * pv # This is there the q-values are calculated
  qv[-1] = min(qv[-1], 1.0)

  # The following ensures that lower ranking q-values are no larger than those
  # that come after
  for i in range(len(pv) - 2, -1, -1):
    qv[i] = min(pi0*m*pv[i] / (i+1.0), qv[i+1])

  qv[p_ordered] = qv.copy() # same as input pv order
  qv = qv.reshape(original_shape) # reshape q-values to orginal shape of pv
  return qv, pi0

# def multipy_est_pi0(pvals, m, s=None, plot=False):
#   kappa = np.arange(0, 0.96, 0.01)
#   pik = [sum(pvals > k) / (m*(1-k)) for k in kappa]
#   cs = UnivariateSpline(kappa, pik, k=3, s=s, ext=0)
#   pi0 = float(cs(1.))
#   print(f'm={len(pvals):<5d} {pi0=:g}')
#   plt.title(f'#tests: {m} {s=} multipy pi0 estimate with spline')
#   plt.plot(kappa, pik, 'o' , markersize=2, label='pi0s')
#   kappa = np.arange(0, 1.01, 0.01)
#   plt.plot(kappa, cs(kappa))
#   plt.plot(1.0, pi0, 'o')
#   plt.show()
#   return pi0

def estimate_pi0(pv, m, pi0=None, plot=False):
  # if the number of hypotheses is small, just set pi0 to 1
  if pi0 is not None:
    return pi0
  elif len(pv) < 100:
    print(f'warning: Too few p-values given ({len(pv)} to estimate pi0. Setting pi0 to 1.0, but this is probably too conservative. You might try using the Benjamini–Hochberg procedure instead. See scipy.stats.false_discovery_control. Also consider the Holm-Bonferroni method: https://en.wikipedia.org/wiki/Holm%E2%80%93Bonferroni_method.', file=sys.stderr)
    return 1.0

  # evaluate pi0 for different lambdas
  pi0s = []
  lam = np.arange(0, 0.95, 0.01)
  lam_len = len(lam)
  counts = np.array([(pv > pvalue).sum() for pvalue in lam])

  for idx in range(lam_len):
    if counts[idx] == 0:
      break
    pi0s.append( counts[idx] / (m * (1 - lam[idx])) )

  pi0s = np.array(pi0s)

  # Estimate pi0 based on the where the pi0s flatten out the most according to
  # the standard deviation. The idea is that the pi0s converge to pi0 as lambda
  # approaches 1. It is probably wise to always inspect the plot to make sure
  # the estimate looks reasonable given the data.
  pi0s_len = len(pi0s)
  if pi0s_len < int(lam_len * 0.75):
    print(f'warning: there were no p-values > ~0.75. Estimate of pi0 is probably overly convervative. It is strongly advised to look at the plot to determine of the pi0 value looks reasonable or provide your own pi0. Also consider another method for estimating the false discovery rate.', file=sys.stderr)
    sd_range = range(pi0s_len - 1)
  else:
    sd_range = range(int(lam_len * 0.6), min(int(lam_len * 0.95), pi0s_len - 1))
  if pi0s_len < lam_len:
    lam = np.delete(lam, range(pi0s_len, lam_len))
  stdev = {}
  means = {}
  min_stdev = np.inf
  min_idx = None
  for idx in sd_range:
    means[idx] = np.mean(pi0s[idx:])
    stdev[idx] = np.std(pi0s[idx:], ddof=1)
    if stdev[idx] < min_stdev:
      min_stdev = stdev[idx]
      min_idx = idx
  pi0 = means[min_idx]

  # print(f'<tr><td>{len(pv)}</td><td>{pi0:g}</td></tr>')
  # print(f'm={len(pv):<5d} {pi0=:g}')
  
  if plot:
    fig, ax1 = plt.subplots()
    plt.title("Estimation of proportion of null hypotheses from p-values");
    ax1.set_xlabel('λ (min p-value for calculating proportion)')
    ax1.set_ylabel(r'pi0s = #p$_j$>λ / m(1-λ) (proportion null)')
    p1 = ax1.plot(lam, pi0s, 'o', markersize=2, label=f'pi0s')
    p2 = ax1.plot(lam[list(means.keys())], list(means.values()), label=f'means')
    p3 = ax1.plot([0, 0.95], [pi0, pi0], label=f'pi0 est={pi0:g}')
    p4 = ax1.plot([lam[min_idx], lam[min_idx]], [pi0, 1], label='min stdev')
    ax2 = ax1.twinx()
    ax2.set_ylabel('std dev of pi0s > λ')
    p5 = ax2.plot(lam[list(stdev.keys())], list(stdev.values()), 'r', label=f'stdevs')
    fig.tight_layout()
    plots = p1+p2+p3+p4+p5
    labels = [plot.get_label() for plot in plots]
    ax1.legend(plots, labels)
    plt.show()


  if pi0 > 1 or pi0 <= 0:
    pi0 = 1.0

  return pi0
