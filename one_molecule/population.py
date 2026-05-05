import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


A=np.loadtxt('corr.txt')
print(A)
t=A[:,0]
t=t/40

B=np.loadtxt('pop.dat')
l=len(B)
t_b=t[0:l]*5

print(A)
print(B)


#######################################################
plt.xlabel('Time (fs)')
plt.ylabel('Population')
###################################################
plt.plot(t,A[:,2],'r-',label=r"S$_1^{'}$",alpha=0.5)
plt.plot(t,A[:,3],'g-',label=r"S$_2^{'}$",alpha=0.5)
plt.plot(t,A[:,4],'b-',label=r"S$_2^{'}$",alpha=0.5)
plt.plot(t,A[:,5],'k-',label=r"S$_3^{'}$",alpha=0.5)
##################################################
plt.plot(t_b,B[:,0],'r--',label=r"S$_1^{'}$(new)")
plt.plot(t_b,B[:,1],'g--',label=r"S$_2^{'}$(new)")
plt.plot(t_b,B[:,2],'b--',label=r"S$_2^{'}$(new)")
plt.plot(t_b,B[:,3],'k--',label=r"S$_3^{'}$(new)")
######################################################
#######################################################
plt.xlim(0,5)
#plt.ylim(-0.1,1.1)
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
plt.savefig('test_pop.png',dpi=400,bbox_inches='tight')
#######################################################

