# qvalue
qvalue function that fixes the issue of misused splines in other versions.

Other python versions of this available on the web do not properly estimate pi0 because it uses splrep with s=0 (no smoothing). The pi0 estimate from these functions is simply the last pi0 value considered and the spline fitting does nothing. e.g.

https://github.com/nfusi/qvalue

https://gist.github.com/ryananeff/0232597b04ec1e5947de2ad8b9292d6e
  
Here, instead of splines, we are calculating the standard deviations of the pi0s > lambda. Since the pi0s tend to flatten out as lambda approaches 1 the standard deviation will tends to be at a minimum at the lambda where the pi0s converge to the null proportion. The pi0 is taken as the mean of the pi0s > lambda where the standard deviation is at a minimum. See https://dxcalc.com/pub/pi0_estimate.png for a graphical illustration.
