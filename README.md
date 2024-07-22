# qvalue
qvalue function that fixes the issue of misused splines in other versions.

Other python versions of this available on the web do not properly estimate pi0 because they use splrep with s=0 (no smoothing). The pi0 estimate from these functions is simply the last pi0 value considered (which is usually a poor estimate) and the spline fitting does nothing. e.g.

https://pypi.org/project/qvalue/

https://github.com/nfusi/qvalue

https://gist.github.com/ryananeff/0232597b04ec1e5947de2ad8b9292d6e
  
Here, instead of splines, we are calculating the standard deviations of the pi0s > lambda. Since the pi0s tend to flatten out as lambda approaches 1 the standard deviation will tend to be at a minimum at the lambda where the pi0s converge to the null proportion. The pi0 is taken as the mean of the pi0s > lambda where the standard deviation is at a minimum. Below is a graphical illustration. The input data here are 28,884 t-test p-values from a differential rna-seq experiment. If the number of input p-values is less then 100 then pi0 is set to 1.0 which is too conservative. For such smaller datasets either determine a reasonable pi0 independently or us another method such as scipy.stats.false_discovery_control. The table below shows the estimates of pi0 when taking random subsamples from the p-values in the plot below. Samples below ~2,000 can give more variable estimates. A dataset with a higher underlying pi0 should be more stable at lower sample sizes.
<table>
  <caption>variance of pi0 estimates as a function of sample size.</caption>
  <th><td># tests</td><td>p0 estimate</td></th>
  <tr><td>m=128</td><td>pi0=0.0786511</td></tr>
  <tr><td>m=256</td><td>pi0=0.307739</td></tr>
  <tr><td>m=512</td><td>pi0=0.388281</td></tr>
  <tr><td>m=1024</td><td>pi0=0.199061</td></tr>
  <tr><td>m=2048</td><td>pi0=0.25974</td></tr>
  <tr><td>m=4096</td><td>pi0=0.255607</td></tr>
  <tr><td>m=8192</td><td>pi0=0.285897</td></tr>
  <tr><td>m=16384</td><td>pi0=0.258258</td></tr>
  <tr><td>m=28884</td><td>pi0=0.255432</td></tr>
</table>

![pi0_estimate](https://github.com/user-attachments/assets/4c54cc9f-8fae-4827-b02c-becf3590e8ca)
