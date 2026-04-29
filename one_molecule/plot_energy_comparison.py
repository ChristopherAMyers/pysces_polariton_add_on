import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


A=np.loadtxt('ex_energy.txt')
print(A)
t=A[:,0]
t=t/40

B=np.loadtxt('polaritonic_energy.dat')
l=len(B)
t_b=t[0:l]


#######################################################
plt.xlabel('Time (fs)')
plt.ylabel('Energy (eV)')
###################################################
plt.plot(t,A[:,1],'r-',label=r"S$_1^{'}$",alpha=0.5)
plt.plot(t,A[:,2],'g-',label=r"S$_2^{'}$",alpha=0.5)
plt.plot(t,A[:,3],'b-',label=r"S$_2^{''}$",alpha=0.5)
plt.plot(t,A[:,4],'k-',label=r"S$_3^{'}$",alpha=0.5)
##################################################
plt.plot(t_b,B[:,0],'r--',label=r"S$_1^{'}$(new)")
plt.plot(t_b,B[:,1],'g--',label=r"S$_2^{'}$(new)")
plt.plot(t_b,B[:,2],'b--',label=r"S$_2^{''}$(new)")
plt.plot(t_b,B[:,3],'k--',label=r"S$_3^{'}$(new)")
######################################################
#######################################################
plt.xlim(0,1)
plt.ylim(2.5,6.8)
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
plt.savefig('test.png',dpi=400,bbox_inches='tight')
#######################################################
