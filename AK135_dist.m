dist=0:6;
tp_z0=[0,19.17,35.03,48.78,62.53,76.27,90.01];
ts_z0=[0,32.14,60.75,85.43,110.10,134.76,159.41];
tpFun=fit(dist',tp_z0','smoothingspline');
tsFun=fit(dist',ts_z0','smoothingspline');
dist2=0:0.001:6;
tpDist2z0=feval(tpFun,dist2');
tsDist2z0=feval(tsFun,dist2');
idx=find(abs(tsMintp-24)==min(abs(tsMintp-24)),1,'first')
dist2(idx)
idx=find(abs(tsMintp-38)==min(abs(tsMintp-38)),1,'first');
dist2(idx)
idx=find(abs(tsMintp-41)==min(abs(tsMintp-41)),1,'first')
dist2(idx)
