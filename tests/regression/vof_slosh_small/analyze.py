import sys, numpy as np
rl, rg, g, h = float(sys.argv[2]), float(sys.argv[3]), 9.81, 0.5
k = 2*np.pi
w = np.sqrt((rl-rg)*g*k/(rl/np.tanh(k*h)+rg/np.tanh(k*h)))
T = 2*np.pi/w
d = np.loadtxt(sys.argv[1])
t = d[:,1]; a = 2*d[:,13]/(1.0*0.078125)  # cos moment column
# amplitude relative to initial; find first and later zero crossings
a0 = a[0]
idx = np.where(np.sign(a[1:])!=np.sign(a[:-1]))[0]
cross = [t[i] - a[i]*(t[i+1]-t[i])/(a[i+1]-a[i]) for i in idx]
print("a0=%.5f theory T=%.5f (half %.5f)"%(a0,T,T/2))
print("zero crossings:", ["%.4f"%c for c in cross])
if len(cross)>=2:
    hp = np.diff(cross)
    print("half periods:", ["%.4f"%x for x in hp], " mean err %.2f%%"%(100*(hp.mean()/(T/2)-1)))
print("steps", len(t), "max |V| col", d[:,15].max())
