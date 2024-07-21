import sys
import scipy as sp
import numpy as np
import matplotlib.pyplot as plt

def qvalue(pv, pi0=None, m=None, verbose=False, plot=False):
  """
  NOTE: Other python versions of this available on the web do not properly
  estimate pi0 because it uses splrep with s=0 (no smoothing). The pi0 estimate
  from these functions is simply the last pi0 value considered and the spline
  fitting does nothing. e.g.
  https://github.com/nfusi/qvalue
  https://gist.github.com/ryananeff/0232597b04ec1e5947de2ad8b9292d6e
  
  Here, instead of splines, we are calculating the standard deviations of the pi0s
  > lambda. Since the pi0s tend to flatten out as lambda approaches 1 the
  standard deviation tends to be at a minimum at the lambda where the
  pi0s converge to the null proportion. The pi0 is taken as the mean of the pi0s
  > lambda where the standard deviation is at a minimum. See
  https://dxcalc.com/pub/pi0_estimate.png for a graphical illustration.

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

  retruns numpy array of q-values with the same shape and order as the input p-values
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

  assert(pi0 >= 0 and pi0 <= 1), f'pi0 is not between 0 and 1: {pi0=:%f}'

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
  return qv

def estimate_pi0(pv, m, pi0=None, plot=False):
  # if the number of hypotheses is small, just set pi0 to 1
  if pi0 is not None:
    return pi0
  elif len(pv) < 100:
    return 1.0

  # evaluate pi0 for different lambdas
  pi0s = []
  lam = np.arange(0, 0.95, 0.01)
  lam_len = len(lam)
  counts = np.array([(pv > pvalue).sum() for pvalue in lam])

  for idx in range(lam_len):
      pi0s.append( counts[idx] / (m * (1 - lam[idx])) )

  pi0s = np.array(pi0s)

  # Estimate pi0 based on the where the pi0s flatten out the most according to
  # the standard deviation. The idea is that the pi0s converge to pi0 as lambda
  # approaches 1. It is probably wise to always inspect the plot to make sure
  # the estimate looks reasonable given the data.
  sd_range = range(int(lam_len * 0.6), int(lam_len * 0.95))
  stdev = {}
  means = {}
  min_stdev = np.inf
  for idx in sd_range:
    means[idx] = np.mean(pi0s[idx:])
    # For some reason sp.std 1.26.2 isn't recognizing the mean argument even
    # though the docs say 1.20. I read that is works with 2.0 though.
    # stdev[idx] = np.std(pi0s[idx:], ddof=1, mean=means[idx])
    stdev[idx] = np.std(pi0s[idx:], ddof=1)
    if stdev[idx] < min_stdev:
      min_stdev = stdev[idx]
      min_idx = idx

  pi0 = means[min_idx]
  # if verbose:
  print(f'qvalue: {pi0=:g}, estimated proportion of null features', file=sys.stderr)

  if plot:
    fig, ax1 = plt.subplots()
    plt.title("Estimation of proportion of null hypotheses from p-values");
    ax1.set_xlabel('λ (min p-value for calculating proportion)')
    ax1.set_ylabel(r'pi0s = #p$_j$>λ / m(1-λ) (proportion null)')
    p1 = ax1.plot(lam, pi0s, 'o', markersize=2, label=f'pi0s')
    p2 = ax1.plot(lam[list(means.keys())], list(means.values()), label=f'means')
    p3 = ax1.plot([0, 0.95], [pi0, pi0], label='pi0 est')
    p4 = ax1.plot([lam[min_idx], lam[min_idx]], [pi0, 1], label='min stdev')
    ax2 = ax1.twinx()
    ax2.set_ylabel('std dev of pi0s > λ')
    p5 = ax2.plot(lam[list(stdev.keys())], list(stdev.values()), 'r', label=f'stdevs')
    fig.tight_layout()
    plots = p1+p2+p3+p4+p5
    labels = [plot.get_label() for plot in plots]
    ax1.legend(plots, labels)
    plt.show()


  if pi0 > 1:
    # if verbose:
    #   print("got pi0 > 1 (%.3f) while estimating q-values, setting it to 1" % pi0)
    pi0 = 1.0

  return pi0
