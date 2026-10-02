import numpy as np,sys
d=np.loadtxt('vof_diag.dat')
print("max -Umin %.3f max|V| %.3f  x(Umax peak) %.2f; Co max %.3f; steps %d; t_end %.2f"%(d[:,14].max(),d[:,15].max(),d[np.argmax(d[:,14]),21],d[:,8].max(),len(d),d[-1,1]))
g=np.loadtxt('vof_gauges.dat'); t=g[:,0]
for T0 in (3,5,7):
    m=(t>T0)&(t<=T0+2)
    if m.sum()>2: print("window",T0,"gauge amps",np.round((g[m,1:].max(0)-g[m,1:].min(0))/2,4))
print("vol drift %.2e bnd %.3e relax %.3e ledger resid %.2e"%(d[-1,3],d[-1,19],d[-1,20],(d[-1,2]-d[0,2]-d[-1,19]-d[-1,20])/d[0,2]))
