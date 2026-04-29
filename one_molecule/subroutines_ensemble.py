import sys
import numpy as np
import subprocess as sp
import random
import pandas
import time
import scipy.integrate as it
#from input_simulation import *
#from input_gamess import option as opt
#from fileIO import SimulationLogger, write_restart, read_restart
# __location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
#__location__ = ''

#from input_simulation import *
natom=10
nel=4
nnuc = 30
ndof = nnuc + nel
init_state=3

# Physical constants and unit conversion factors
pi       = np.pi
hplanck  = 6.62607015*10**-34   # Planck's constant in SI
hbar     = hplanck/(2*pi)       # reduced Planck's constant in SI
clight   = 2.99792458*10**8     # speed of light in SI
kb       = 1.380649*10**-23     # Boltzmann constant
eh2j     = 4.359744650*10**-18  # Hartree to Joule
amu2au   = 1.822888486*10**3    # atomic mass unit to atomic unit
autime2s = hbar/eh2j            # atomic unit time to second
au2fs    = autime2s*10**15      # atomic unit time to second
ang2bohr = 1.8897259886         # angstroms to bohr
k2autmp  = kb/eh2j              # Kelvin to atomic unit temperature
#beta     = 1.0/(temp * k2autmp) # inverse temperature in atomic unit



# Get atom labels
def get_atom_label():
    atoms = []
    if(mol_input_format == "gamess"):
      f = open(os.path.join(__location__,'geo_gamess'), 'r')
      f.readline()
    elif(mol_input_format == "terachem"):
      f = open(os.path.join(__location__,fname_tc_xyz), 'r')
      f.readline()
      f.readline()
    for i in range(natom):
        x = f.readline().split()
        atoms.append(x[0])
    return(atoms)



########################################################
#####Nuclear phase space variables######################
######################################################
def sample_nuclear_nm(qcenter):
    Q = qcenter
    P = pN0
    return(Q, P)






'''LSC-IVR with Wigner population estimator'''
def sample_wignerLSC():
    coord = np.zeros((2, ndof-6))

    # Determine the sampling radius of initially occupied electronic state
    from scipy.optimize import fsolve
    from functools import partial
    def eqn(F, r):
        return 2 ** (F + 1) * (r - 0.5) * np.exp(-(r + 0.5 * (F - 1))) - 1

    root = fsolve(partial(eqn, nel), 2)
    r = root[0] ** 0.5

    # Electronic phase space variables
    for i in range(nel):
        theta = random.random()
        if i == init_state-1:
            x = r * np.cos(2.0*pi*theta)
            p = r * np.sin(2.0*pi*theta)
        else:
            x = np.sqrt(1.0/2.0) * np.cos(2.0*pi*theta)
            p = np.sqrt(1.0/2.0) * np.sin(2.0*pi*theta)

        coord[0, i] = x
        coord[1, i] = p

# Nuclear phase space variables
#    for i in range(nnuc-6):
#       coord[0, i+nel], coord[1, i+nel] = sample_nuclear(qN0[i], frq[i+6])

#    # Nuclear phase space variables
#    for i in range(nnuc-6):
#        # position
#        coord[0,i+nel] = np.random.normal(loc=qN0[i], scale=np.sqrt(1.0/(2.0*frq[i+6]*np.tanh(beta*frq[i+6]/2.0))))
#        # momentum
#        coord[1,i+nel] = np.random.normal(loc=pN0, scale=np.sqrt(frq[i+6]/(2.0*np.tanh(beta*frq[i+6]/2.0))))

    return(coord)





#####################################################################
### Compute equations of motion (mapping variables derivatives)   ###
### of adiabatic MM-ST Hamiltonian with the symmetrized potential ###
#####################################################################
def get_derivatives(au_mas, q, p, nac, grad, elecE):
    der = np.zeros((2, ndof))

    # Derivatives of elctronic mapping variables
    for i in range(nel):
        xdpm   = 0 # NAC part of the derivative, x*d*p/m
        pdpm   = 0 # NAC part of the derivative, p*d*p/m
        sum_DE = 0 # Traceless Hamiltonian, sum(Ei - Ej)
        for j in range(nel):
            if j != i:
                xdpm   += q[j] * np.matmul(nac[j,i,:], p[nel:]/au_mas)    # NAC part of position derivative
                pdpm   += p[j] * np.matmul(nac[j,i,:], p[nel:]/au_mas)    # NAC part of momenta derivative
                sum_DE += elecE[i] - elecE[j]
        # postions
        der[0, i] =  (1.0/nel) * p[i] * sum_DE + xdpm     # total position derivative
        # momenta
        der[1, i] = -(1.0/nel) * q[i] * sum_DE + pdpm     # total momenta derivative
   # momenta
        diag_part = -(1.0/nel) * q[i] * sum_DE
        nac_part  = pdpm
        tmp_pdpm = p[j] * np.matmul(nac[j,i,:], p[nel:]/au_mas)
        print("i,j =", i, j, "pdpm_contrib =", tmp_pdpm)
#der[1, i] = diag_part + nac_part

#        print(i, "dp_diag =", diag_part, "dp_nac =", nac_part, "dp_total =", der[1, i])
    
    # Derivatives of nuclear mapping variables
    for n in range(nnuc):
        # positions
        der[0, nel+n] = p[nel+n]/au_mas[n]
        # momenta
        sum_dEdR   = 0
        p2x2_DdEdR = 0 # (pi^2 - pj^2 + qi^2 - qj^2) * (dEi/dR - dEj/dR)
        ppxx_DEnac = 0 # (pi * pj + qi * qj) * (Ej - Ei) * dij
        for i in range(nel):
            sum_dEdR  += grad[i,n]
            j = i + 1
            while j < nel:
                p2x2_DdEdR += (p[i]**2 - p[j]**2 + q[i]**2 - q[j]**2) * (grad[i,n] - grad[j,n])
                ppxx_DEnac += (p[i]*p[j] + q[i]*q[j]) * (elecE[j] - elecE[i]) * nac[i,j,n]
                j += 1
        der[1, nel+n] = -(1.0/nel) * sum_dEdR - (0.5/nel) * p2x2_DdEdR - ppxx_DEnac
    # momenta


    return(der)



# Compute total energies to verify that it is conserved
def get_energy(au_mas, q, p, elecE):
    # Nuclear part (sum of P**2/M)
    p2m_sum = 0
    for n in range(nnuc):
        p2m_temp = 0
        p2m_temp = p[nel+n]**2/au_mas[n]
        p2m_sum += p2m_temp

    # Electronic part ((pi2 - p2j + qi2 - qj2) * (Ei - Ej))
    p2x2_DE = 0
    for i in range(nel):
        j = i + 1
        while j < nel:
            p2x2_DE += (p[i]**2 - p[j]**2 + q[i]**2 - q[j]**2) * (elecE[i] - elecE[j])
            j += 1

    # Total energy at updated t
    energy = 0.5*p2m_sum + (1.0/nel)*sum(elecE) + (0.5/nel)*p2x2_DE
    return(energy)


# Simple scipy based integrator
def scipy_rk4(elecE, grad, nac, yvar, dt, au_mas):
    def get_deriv(t, y0):
        der = get_derivatives(au_mas, y0[:ndof], y0[ndof:], nac, grad, elecE)
        der = der.flatten()
        return(der)
    result = it.solve_ivp(get_deriv, (0,dt), yvar, method='RK45', max_step=dt, t_eval=[dt], rtol=1e-10, atol=1e-10)
    return(result.y.flatten())




# Actual RK4 integrator
def integrate_rk4(elecE, grad, nac, dt, yvar, au_mas):
#   au_mas = np.diag(amu_mat) * amu2au
   h  = dt
   y0 = yvar.copy()
   y1 = np.zeros(2*ndof)
   y2 = np.zeros(2*ndof)
   y3 = np.zeros(2*ndof)

   ### 4th-order Runge-Kutta routine ###
   # Get derivatives (k1) at y0
   k1 = get_derivatives(au_mas, y0[:ndof], y0[ndof:], nac, grad, elecE)
   k1 = k1.flatten()
   for i in range(2*ndof):
        y1[i] = y0[i] + 0.5*h*k1[i]

   # Get derivatives (k2) at y1
   k2 = get_derivatives(au_mas, y1[:ndof], y1[ndof:], nac, grad, elecE)
   k2 = k2.flatten()
   # Take an intermediate step to y2
   for i in range(2*ndof):
        y2[i] = y0[i] + 0.5*h*k2[i] # position

   # Get derivatives (k3) at y2
   k3 = get_derivatives(au_mas, y2[:ndof], y2[ndof:], nac, grad, elecE)
   k3 = k3.flatten()
   # Take an intermediate step to y3
   for i in range(2*ndof):
        y3[i] = y0[i] + h*k3[i] # position

   # Get derivatives (k4) at y3
   k4 = get_derivatives(au_mas, y3[:ndof], y3[ndof:], nac, grad, elecE)
   k4 = k4.flatten()

   # Compute the coordinates at t=x+H
   result = np.zeros(2*ndof)
   for i in range(2*ndof):
      result[i] = y0[i] + h*(k1[i] + 2*k2[i] + 2*k3[i] + k4[i])/6.0

   return(result)











