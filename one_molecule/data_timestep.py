import numpy as np
import os


# Check NAC sign
def correct_nac_sign(nac, nac_old):
    for i in range(0,n_elec):
        for j in range(0,n_elec):
            if (i != j):
                dotval = np.dot(nac[i,j,:], nac_old[i,j,:])
                sign = np.sign(dotval)
                if (sign == 0):
                    sign = 1.0
                nac[i,j,:] = sign*nac[i,j,:]
    return nac



# Main path
path_data=os.getcwd()
print ("Main path = ", path_data)


# Test
n_m=1                 #np.loadtxt(n_m.dat)
n_m=int(n_m)
n_state=3
n_elec=n_state+1
n_nuc=10


# initialize
energy=np.zeros([n_m,n_elec])
gradient=np.zeros([n_m,n_elec,n_nuc,3])
na_coupling=np.zeros([n_m,n_elec,n_elec,n_nuc,3])
mu=np.zeros([n_m,n_elec,n_elec,3])
tdp=np.zeros([n_m,n_elec,n_elec,n_nuc,3,3])




# load all data
# energies
for i in range(0, n_m):
    path_traj = os.path.join(path_data.strip(),f'traj_{i+1}')
    print('Trajectory (energies) : ',path_traj)
    energy[i,0]=0.0
    for j in range(1,n_elec):
        E0=np.loadtxt(os.path.join(path_traj.strip(),'gs.dat'))
        E=np.loadtxt(os.path.join(path_traj.strip(),'ex.dat'))
        energy[i,j]=E[j-1] - E0  # /27.21140795
# gradients and nacs
for i in range(0, n_m):
    path_traj = os.path.join(path_data.strip(),f'traj_{i+1}')
    print('Trajectory (grads and nacs) : ',path_traj)
#############################################################################    
    for j in range(0,n_elec):
        path = os.path.join(path_traj.strip(),f'grad_{j}.dat')
        grad = np.loadtxt(path)
        print(j)
        print(grad)
        gradient[i,j,:,:]=grad
#############################################################################        
        for k in range(j+1,n_elec):
            path = os.path.join(path_traj.strip(),f'nac_{j}{k}.dat')
            nac = np.loadtxt(path)
            na_coupling[i,j,k,:,:]=(+1)*nac
            na_coupling[i,k,j,:,:]=(-1)*nac
            print(j,k)
            print(nac)
#############################################################################
# dipole moments
for i in range(0, n_m):
    path_traj = os.path.join(path_data.strip(),f'traj_{i+1}')
    print('Trajectory (dipole moment) : ',path_traj)
    path = os.path.join(path_traj.strip(),'u.dat')
    u=np.loadtxt(path)
    for j in range(0,1):
        for k in range(j+1,n_elec):
            mu[i,j,k,:] = u[k-1,:] 
            mu[i,k,j,:] = u[k-1,:]
            print(i,j,k)
            print(u[k-1,:])
#############################################################################
# transition dipole moment derivatives
#############################################################################
for i in range(0, n_m):
    path_traj = os.path.join(path_data.strip(),f'traj_{i+1}')
    print('Trajectory (transition dipole derivatives) : ',path_traj)
#######################################
    for j in range(0,1):
        for k in range(j+1,n_elec):
########################################
            path_x = os.path.join(path_traj.strip(),f'tdp_0{k}_x.dat')
            path_y = os.path.join(path_traj.strip(),f'tdp_0{k}_y.dat')
            path_z = os.path.join(path_traj.strip(),f'tdp_0{k}_z.dat')
##########################################
            tdp_x = np.loadtxt(path_x)
            tdp_y = np.loadtxt(path_y)
            tdp_z = np.loadtxt(path_z)
##########################################
            tdp[i,j,k,:,0,:] = tdp_x
            tdp[i,j,k,:,1,:] = tdp_y
            tdp[i,j,k,:,2,:] = tdp_z
###########################################
            tdp[i,k,j,:,0,:] = tdp_x
            tdp[i,k,j,:,1,:] = tdp_y
            tdp[i,k,j,:,2,:] = tdp_z
###########################################
            print(i,j,k)
            print(tdp_x)
            print(tdp_y)
            print(tdp_z)
#############################################################################
################################################################################



#################################################################################
# Do sign check with NAC from last step
nac_flat = na_coupling.reshape(n_m,n_elec,n_elec,n_nuc*3)

if os.path.exists("nac_hist.npy"):
    nac_old = np.load("nac_hist.npy")

    for i in range(0,n_m):
        nac_flat[i,:,:,:] = correct_nac_sign(nac_flat[i,:,:,:], nac_old[i,:,:,:])

nac_old = nac_flat.copy()

na_coupling = nac_flat.reshape(n_m,n_elec,n_elec,n_nuc,3)

np.save("nac_hist.npy", nac_old)
####################################################################################



# save data
print('Saving Data')
np.save("tdp.npy", tdp)
np.save("gradient.npy", gradient)
np.save("nac.npy", na_coupling)
np.save("mu.npy", mu)
np.save('energy.npy', energy)
print('Done')
