import sys
import numpy as np
import matplotlib.pyplot as plt

def qvalue(pv, pi0=None, m=None, plot=False, verbose=False):
  """Estimates q-values from p-values based on Storey and Tibshirani, 2003
  Args
  =====
  pv:      array of p-values
  pi0:     If None, it's estimated based on principles discussed in Storey and Tibshirani, 2003, but see below.
  m:       Number of tests. If not specified m = pv.size
  plot:    If True and pi0=None, display a plot showing how pi0 was estimated.
  verbose: Print verbose messages? Not currently used. (default False)

  returns qv, pi0
    qv:  a numpy array of q-values with the same shape and order as the input p-values
    pi0: the pi0 estimate

  This was written because other python versions of this function were not
  properly estimating pi0. See
  https://github.com/Old-Green-Man/qvalue/blob/main/README.md for more
  information.

  It is probably a good idea to call this with plot=True on your p-values, at
  least once, to see what the distribution of pi0s looks like. If for some
  reason it seems off you can try providing your own pi0 with the pi0 function
  option. A higher pi0 will give more conservative q-value estimates. The plot
  should give you an idea if you can make a reasonable estimate of pi0 or
  not. Some sets of p-values can give noisy plots where it's hard to know what
  the pi0 should be. In such cases you could try a very conservative (higher)
  pi0 or try another FDR method such as scipy.stats.false_discovery_control.

  """

  pv = np.asarray(pv)
  assert(pv.min() >= 0 and pv.max() <= 1), "p-values should be between 0 and 1"

  original_shape = pv.shape
  pv = pv.ravel() # flatten

  pv_len = len(pv)
  if m is None:
    m = pv_len
  else: # user specified
    m *= 1.0

  pi0 = estimate_pi0(pv, m, pi0, plot=plot)
  # for pwr in range(-10, 11):
  #   Note: s=None is used in multipy
  #   pi0 = multipy_est_pi0(pv, m, s=2**pwr, plot=plot)
  # pi0 = multipy_est_pi0(pv, m, plot=plot)
  

  p_ordered = np.argsort(pv)
  pv = pv[p_ordered]
  qv = np.zeros_like(pv)
  qv[-1] = pi0 * m * pv[-1] / pv_len

  # The following min(x, y) ensures that lower ranking q-values are no larger
  # than those that come after.

  # In the equation below for the qvalue, (pi0 * m * pv[i]) / (i+1.0), the
  # numerator is the expected number of false positives, since pi0*m is the
  # estimated number of truly null hypotheses and pv[i] is the probability of a
  # truly null feature being called significant (being below the threshold
  # t). The denominator is the number of features called significant.  The
  # q-value for a feature then is the minimum FDR that can be attained when
  # calling that feature significant.
  for i in range(pv_len - 2, -1, -1):
    qv[i] = min(pi0 * m * pv[i] / (i+1.0), qv[i+1])

  qv[p_ordered] = qv.copy() # same as input pv order
  qv = qv.reshape(original_shape) # reshape q-values to orginal shape of pv
  return qv, pi0

def estimate_pi0(pv, m, pi0=None, plot=False):
  # if the number of hypotheses is small, just set pi0 to 1
  pv_len = len(pv)
  if pi0 is not None:
    return pi0
  elif pv_len < 100:
    print(f'warning: Too few p-values given ({pv_len} to estimate pi0. Setting pi0 to 1.0, but this is probably too conservative. You might try using the Benjamini–Hochberg procedure instead. See scipy.stats.false_discovery_control. Also consider the Holm-Bonferroni method: https://en.wikipedia.org/wiki/Holm%E2%80%93Bonferroni_method.', file=sys.stderr)
    return 1.0

  # evaluate pi0 for different lambdas
  pi0s = []
  # This range should probably be adjusted so that # p-values > lambda is above some minimum
  lam = np.arange(0, 0.95, 0.01)
  lam_len = len(lam)
  counts = np.array([(pv > pvalue).sum() for pvalue in lam])
  for idx in range(lam_len):
    if counts[idx] == 0: # stop when we have no pvalues larger than this lambda
      break
    pi0s.append( counts[idx] / (m * (1 - lam[idx])) )

  pi0s = np.array(pi0s)

  # Estimate pi0 based on the where the pi0s flatten out the most according to
  # the standard deviation. The idea is that the pi0s converge to pi0 as lambda
  # approaches 1. It is probably wise to always inspect the plot to make sure
  # the estimate looks reasonable given the data.
  pi0s_len = len(pi0s)

  # The commented out section below was where I was trying to start at higher
  # lambda values so I wouldn't have to calculate as many standard
  # deviations. However, it's probably safer to just consider all the
  # lambdas. Samples that have very few non-null tests (pi0 close to 1) will
  # have better fits toward lambda = 0.
  
  # if pi0s_len < int(lam_len * 0.75):
  #   print(f'warning: there were no p-values > ~0.75. Estimate of pi0 is probably overly convervative. It is strongly advised to look at the plot to determine of the pi0 value looks reasonable or provide your own pi0. Also consider another method for estimating the false discovery rate.', file=sys.stderr)
  #   sd_range = range(pi0s_len - 1)
  # else:
  #   sd_range = range(int(lam_len * 0.6), min(int(lam_len * 0.95), pi0s_len - 1))
  
  sd_range = range(min(int(lam_len * 0.90), pi0s_len - 5))
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

  # print(f'<tr><td>{pv_len}</td><td>{pi0:g}</td></tr>')
  # print(f'm={pv_len:<5d} {pi0=:g}')
  
  if plot:
    fig, ax1 = plt.subplots()
    plt.title("Estimation of proportion of null hypotheses from p-values");
    ax1.set_xlabel('λ (min p-value for calculating proportion)')
    ax1.set_ylabel(r'pi0s = #p$_j$>λ / m(1-λ) (proportion null)')
    p1 = ax1.plot(lam, pi0s, 'o', markersize=2, label=f'pi0s')
    p2 = ax1.plot(lam[list(means.keys())], list(means.values()), label=f'means')
    # p3 = ax1.plot([0, 0.95], [pi0, pi0], label=f'pi0 est={counts[min_idx]}/({m:.0f}*(1-{lam[min_idx]}))={pi0:g}')
    p3 = ax1.plot([0, 0.95], [pi0, pi0], label=f'pi0 est={pi0:g}')
    p4 = ax1.plot([lam[min_idx], lam[min_idx]], [pi0, 1], color='c', label='min stdev')
    ax2 = ax1.twinx()
    ax2.set_ylabel('std dev of pi0s > λ')
    p5 = ax2.plot(lam[list(stdev.keys())], list(stdev.values()), 'r', label=f'stdevs')
    fig.tight_layout()
    plots = p1+p2+p3+p4+p5
    labels = [plot.get_label() for plot in plots]
    ax1.legend(plots, labels)
    plt.show()


  if pi0 > 1: # or pi0 <= 0:
    pi0 = 1.0

  return pi0

# For testing the spline method of estimating pi0
def multipy_est_pi0(pvals, m, s=None, plot=False):
  from scipy.interpolate import UnivariateSpline
  kappa = np.arange(0, 0.96, 0.01)
  pik = [sum(pvals > k) / (m*(1-k)) for k in kappa]
  cs = UnivariateSpline(kappa, pik, k=3, s=s, ext=0)
  pi0 = float(cs(1.))
  print(f'm={len(pvals):<5d} {pi0=:g}')
  plt.title(f'# tests: {m:.0f} {s=} multipy {pi0=:g} estimate with spline')
  plt.plot(kappa, pik, 'o' , markersize=2, label='pi0s')
  kappa = np.arange(0, 1.01, 0.01)
  plt.plot(kappa, cs(kappa))
  plt.plot(1.0, pi0, 'o')
  plt.show()
  return pi0
