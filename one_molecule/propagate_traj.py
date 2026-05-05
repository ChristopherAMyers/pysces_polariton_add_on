import numpy as np
from subroutines_ensemble import sample_wignerLSC, scipy_rk4, get_energy, ndof, nel, nnuc, integrate_rk4, get_derivatives
#from main_test import build_inputs
from coupled_ensemble import CoupledMolecule
import scipy.integrate as it
from ase.io import read
import numpy as np
import subprocess
import os


# read coords and momenta using ASE
def read_xyz(fname):
    atoms = read(fname)
    coords = atoms.get_positions()   # shape (n_atoms, 3)
    labels = atoms.get_chemical_symbols()
    return coords.reshape(-1), labels


# Reshape coordinates to column vector
def reshape(A):
    m,n = np.shape(A)
    coord = np.zeros([m*n])
    index = 0
    for i in range(0,m):
       for j in range(0,n):
           coord[index] = A[i,j]
           index = index + 1
    return(coord) 


# Write updated coordinates
def write_xyz(fname, vec, labels, comment):
    n_atoms = len(labels)
    xyz = vec.reshape(n_atoms, 3)
    with open(fname, "w") as f:
        f.write(f"{n_atoms}\n")
        f.write(f"{comment}\n")
        for i in range(n_atoms):
            f.write(
                f"{labels[i].strip()}  {xyz[i,0]:16.10f}  {xyz[i,1]:16.10f}  {xyz[i,2]:16.10f}\n"
            )


# Compute electronic state population from wigner distribution
def get_pop_wigner(q_e, p_e):
    pop = np.zeros(len(q_e))
    distribution = 2**(len(q_e)+1) * np.exp(-np.dot(q_e, q_e) - np.dot(p_e, p_e))
    for i in range(len(q_e)):
        pop[i] = (q_e[i]**2 + p_e[i]**2 - 0.5)*distribution
    return pop


# Save population and energies to file
def append_vector(fname, vec):
    with open(fname, "a") as f:
        for x in vec:
            f.write(f"{x:16.10f}  ")
        f.write("\n")



# Build inputs for PYSCES
def build_inputs(ref_eig_vecs):
# Physical parameters
    omega_c = 5.55/27.21140795
    gc = 0.0673174/27.21140795
    n_m = 1
    n_elec = 4
    n_nuc = 10
    field_dir = [1.0, 0.0, 0.0]
# States considered (singly excited manifold)
    mol_grads = [[1,2,3]]
# Load Terachem calculation results    
    mol_energies = np.load("energy.npy")
    mol_gradients_raw = np.load("gradient.npy")
    mol_NACs_raw = np.load("nac.npy")
    mol_dipole_matrix = np.load("mu.npy")
    mol_dipole_matrix_gradient_raw = np.load("tdp.npy")
# Reshape
    mol_gradients = mol_gradients_raw.reshape(n_m, n_elec, n_nuc * 3)
    mol_NACs = mol_NACs_raw.reshape(n_m, n_elec, n_elec, n_nuc * 3)
#    mol_dipole_matrix_gradient = np.zeros((n_m, n_elec, n_elec, n_nuc * 3, 3))
    mol_dipole_matrix_gradient = mol_dipole_matrix_gradient_raw.reshape(n_m,n_elec,n_elec,n_nuc*3,3)
# Set up coupled molecule object
    system = CoupledMolecule(
        omega_c=omega_c,
        mol_grads=mol_grads,
        n_elec=n_elec,
        n_nuc=n_nuc,
        n_m=n_m,
        field_dir=field_dir,
        gc=gc,
    )
# Compute polaritonic quantities
    system.compute_all(
        mol_energies=mol_energies,
        mol_gradients=mol_gradients,
        mol_NACs=mol_NACs,
        mol_dipole_matrix=mol_dipole_matrix,
        mol_dipole_matrix_gradient=mol_dipole_matrix_gradient,
        ref_eig_vecs=ref_eig_vecs,
    )

    # system.eigen_val_gradients is gradient of the relative polariton energy
# Add back grad(E0) so nuclear force uses the absolute polariton surface
    ground_gradient = np.load("ground_gradient.npy")
    ground_gradient = ground_gradient.reshape(n_m * n_nuc * 3)

    pol_grad = system.eigen_val_gradients.copy()
    pol_grad = pol_grad + ground_gradient[None, :]

# Return polaritonic stuff    
    return system.eigen_vals, pol_grad, system.NACs, system.eigen_vecs  #system.eigen_val_gradients



# Main driver
def main():
# Time step size    
    dt = 4.0
    nel=4
# Atomic masses
    au_mas = np.loadtxt("mass_new")
    au_mas = au_mas*1.822888486*10**3          #*1822.88848
#    amu_mat = np.diag(au_mas)
#Load nuclear positions and momenta
    q_nuc, labels = read_xyz("q1.xyz")
    p_nuc, _      = read_xyz("p1.xyz") 
    q_nuc = q_nuc*1.8897259886             #1.88972
#    p_nuc = p_nuc*1.88972
# Storing all phase space variables together
    line = subprocess.check_output("tail -n 1 px.dat", shell=True, text=True).strip()
    vals = np.array(line.split(), dtype=float)
    # vals = [time, p0, q0, p1, q1, ...]
    mapvals = vals[1:]   # skip time
    p_el = mapvals[0::2].copy()
    q_el = mapvals[1::2].copy()
    q = np.zeros(ndof)
    p = np.zeros(ndof)
    q[0:nel] = q_el
    p[0:nel] = p_el
    q[nel:] = q_nuc
    p[nel:] = p_nuc
# Total vector to be propagated just as PYSCES wants it
    y = np.concatenate((q, p))
# Build PYSCES inputs from terahcem results and hamiltonian constructtion
#    elecE, grad, nac, eigvecs = build_inputs(ref_eig_vecs=None)
    if os.path.exists("eigvecs_old.npy"):
        ref_eig_vecs = np.load("eigvecs_old.npy")
    else:
        ref_eig_vecs = None
    elecE, grad, nac, eigvecs = build_inputs(ref_eig_vecs=ref_eig_vecs)
    np.save("eigvecs_old.npy", eigvecs)

# Verify matrix and vector shapes    
    print("q shape   =", q.shape)
    print("p shape   =", p.shape)
    print("grad shape=", grad.shape)
    print("nac shape =", nac.shape)
#    E0 = get_energy(au_mas, q, p, elecE)
#    print("Initial energy =", E0)
# Compute and store electronic state populations
    q_old = y[0:ndof]
    p_old = y[ndof:]
###################################
    q_e_old = q_old[0:nel]
    p_e_old = p_old[0:nel]
####################################
    pop = get_pop_wigner(q_e_old, p_e_old)
    elecE_ev = elecE * 27.21140795
    elecE_ev = np.sort(elecE_ev)
###########################################
    append_vector("pop.dat", pop)
    append_vector("polaritonic_energy.dat", elecE_ev)
# Propagate
    der = get_derivatives(au_mas, q, p, nac, grad, elecE)

    print("elecE eV =", elecE*27.21140795)
    print("dq_e/dt =", der[0,0:nel])
    print("dp_e/dt =", der[1,0:nel])

    ref_eig_vecs = eigvecs.copy()    
#    y = scipy_rk4(elecE, grad, nac, y, dt, au_mas)
    y = integrate_rk4(elecE, grad, nac, dt, y, au_mas)
    np.savetxt("traj_status", y)
# Save updated coordinates
    q_new = y[0:ndof]
    p_new = y[ndof:]
########################################    
    q_cart_new = q_new[nel:]
    p_cart_new = p_new[nel:]
#######################################
    n_atom = len(labels)
    n_cart = 3*n_atom
###############################################
    q_mol1 = q_cart_new[0:n_cart]/1.8897259886         #* 0.52917
#    q_mol2 = q_cart_new[n_cart:2*n_cart]
################################################
    p_mol1 = p_cart_new[0:n_cart]   # *0.52917
#    p_mol2 = p_cart_new[n_cart:2*n_cart]
################################################
    write_xyz("q1.xyz", q_mol1, labels, "positions molecule 1")
#    write_xyz("q_mol2.xyz", q_mol2, labels, "positions molecule 2")
#################################################
    write_xyz("p1.xyz", p_mol1, labels, "momenta molecule 1")
#    write_xyz("p_mol2.xyz", p_mol2, labels, "momenta molecule 2")
#Update electronic populations in file
    q_e_new = q_new[0:nel]
    p_e_new = p_new[0:nel]
    time_now = vals[0] + dt
    with open("px.dat", "a") as f:
        f.write(f"{time_now:16.10f}  ")
        for i in range(nel):
            f.write(f"{p_e_new[i]:16.10f}  {q_e_new[i]:16.10f}  ")
        f.write("\n")


if __name__ == "__main__":
    print ("Number of electronic states = ",nel)
    main()
#    print(nel)
