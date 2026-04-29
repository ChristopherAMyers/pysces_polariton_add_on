# -*- coding: utf-8 -*-
"""
Created on Sat Mar  7 20:02:57 2026

@author: Sourav
"""

import numpy as np
from copy import deepcopy
from adiabatic_ensemble import AdiabaticStates

class CoupledMolecule(AdiabaticStates):

    def __init__(self,
                 omega_c: float,
                 mol_grads: list[list[int]],
                 n_elec: int,
                 n_nuc: int,
                 n_m: int,
                 field_dir: list[float],
                 gc: float,
                 ):
        # active states for each molecule
        self.mol_grad_indices = tuple(tuple(states) for states in mol_grads)

        self._omega_c = omega_c
        self._gc = gc
        self._n_elec = n_elec
        self._n_mol = n_m
        self._n_nuc_mol = n_nuc

        # normalize cavity field direction
        self._field_dir = np.asarray(field_dir, dtype=float)
        self._field_dir = self._field_dir / np.linalg.norm(self._field_dir)

        # --------------------------------------------------
        # basis: |e_a^(m),0> and |g,1>
        # --------------------------------------------------

        state_pairs = []

        for m in range(n_m):
            for a in self.mol_grad_indices[m]:
                if a == 0:
                    continue
                state_pairs.append((m, a, 0))

        state_pairs.append((-1, 0, 1))   # photonic state

        self._state_pairs = tuple(state_pairs)
        self._n_dim = len(self._state_pairs)

        super().__init__(self._n_dim, n_m * n_nuc)

        # --------------------------------------------------
        # n_m monomer properties
        # --------------------------------------------------
        self.mol_energies = np.zeros((self._n_mol, self._n_elec))
        self.mol_gradients = np.zeros((self._n_mol, self._n_elec, n_nuc * 3))
        self.mol_NACs = np.zeros((self._n_mol, self._n_elec, self._n_elec, n_nuc * 3))
        self.mol_dipole_matrix = np.zeros((self._n_mol, self._n_elec, self._n_elec, 3))
        self.mol_dipole_matrix_gradient = np.zeros(
            (self._n_mol, self._n_elec, self._n_elec, n_nuc * 3, 3)
        )

    def copy(self):
        return deepcopy(self)

    def compute_all(self,
                    mol_energies,
                    mol_gradients,
                    mol_NACs,
                    mol_dipole_matrix,
                    mol_dipole_matrix_gradient,
                    ref_eig_vecs):
        self.mol_energies = np.copy(mol_energies)
        self.mol_gradients = np.copy(mol_gradients)
        self.mol_NACs = np.copy(mol_NACs)
        self.mol_dipole_matrix = np.copy(mol_dipole_matrix)
        self.mol_dipole_matrix_gradient = np.copy(mol_dipole_matrix_gradient)

        self.set_hamiltonian(self.mol_energies, self.mol_dipole_matrix)
        self.set_hamiltonian_gradient(
            self.mol_gradients,
            self.mol_dipole_matrix,
            self.mol_dipole_matrix_gradient
        )

        self.diagonalize_H(ref_eig_vecs)

        self.calc_NA_coupling(self.mol_NACs)
        self.calc_eigen_value_gradient()

    def calc_NA_coupling(self, mol_basis_NACs):
        basis_NACs = self.get_basis_NACs(mol_basis_NACs)
        return super().calc_NA_coupling(basis_NACs)

    # mu.E
    def _calc_mu_dot_field(self, dipole_matrix):
        n_elec = self._n_elec
        mu_dot_field = np.zeros((n_elec, n_elec))
        for a in range(n_elec):
            for b in range(n_elec):
                mu_dot_field[a, b] = np.dot(dipole_matrix[a, b], self._field_dir)
        return mu_dot_field

    # dipole gradient in direction of cavity field
    def _calc_grad_mu_dot_field(self, dipole_matrix, dipole_matrix_gradient):
        n_elec = self._n_elec
        n_cart = self._n_nuc_mol * 3

        grad_mu_dot_field = np.zeros((n_elec, n_elec, n_cart))
        for a in range(n_elec):
            for b in range(n_elec):
                for k in range(n_cart):
                    grad_mu_dot_field[a, b, k] = np.dot(
                        dipole_matrix_gradient[a, b, k],
                        self._field_dir
                        )
        return grad_mu_dot_field

    def set_hamiltonian(self, mol_energies, mol_dipole_matrix):
        '''
        Evaluate the ensemble Jaynes-Cummings Hamiltonian in the localized
        single-excitation basis:

        |e_a^(m),0>
        |g,1>
        '''
        H_t = np.zeros((self._n_dim, self._n_dim))

        ph = self._state_pairs.index((-1, 0, 1))

        H_t[ph, ph] = self._omega_c

        for i in range(self._n_dim):
            if i == ph:
                continue

            m, a, n = self._state_pairs[i]
            H_t[i, i] = mol_energies[m, a]

        for i in range(self._n_dim):
            if i == ph:
                continue

            m, a, n = self._state_pairs[i]

            mu_dot_field = self._calc_mu_dot_field(mol_dipole_matrix[m])
            gma = self._gc * mu_dot_field[0, a]

            H_t[ph, i] = gma
            H_t[i, ph] = gma

        self.mol_energies = np.copy(mol_energies)
        self.mol_dipole_matrix = np.copy(mol_dipole_matrix)
        self.hamiltonian = H_t
        return H_t

    def set_hamiltonian_gradient(self, mol_gradients, mol_dipole_matrix, mol_dipole_matrix_gradient):
        '''
        Gradient of the ensemble Jaynes-Cummings Hamiltonian in the localized
        single-excitation basis:

        |e_a^(m),0>
        |g,1>
        '''
        n_cart_tot = self.n_nuclei * 3
        n_cart_mol = self._n_nuc_mol * 3

        dH = np.zeros((self._n_dim, self._n_dim, n_cart_tot))

        ph = self._state_pairs.index((-1, 0, 1))

        grad_mu_proj = [self._calc_grad_mu_dot_field(mol_dipole_matrix[m], mol_dipole_matrix_gradient[m])
                    for m in range(self._n_mol)]

        dH[ph, ph, :] = 0.0

        for i in range(self._n_dim):
            if i == ph:
                continue

            m, a, n = self._state_pairs[i]

            i0 = m * n_cart_mol
            i1 = (m + 1) * n_cart_mol

            dH[i, i, i0:i1] = mol_gradients[m, a, :]

        for i in range(self._n_dim):
            if i == ph:
                continue

            m, a, n = self._state_pairs[i]

            i0 = m * n_cart_mol
            i1 = (m + 1) * n_cart_mol

            dg = self._gc * grad_mu_proj[m][0, a, :]

            dH[ph, i, i0:i1] = dg
            dH[i, ph, i0:i1] = dg

        self.mol_gradients = np.copy(mol_gradients)
        self.mol_dipole_matrix = np.copy(mol_dipole_matrix)
        self.mol_dipole_matrix_gradient = np.copy(mol_dipole_matrix_gradient)
        self.dH = dH
        return dH

    def get_basis_NACs(self, mol_basis_NACs):
        '''
        Build diabatic/monomer basis NACs in the localized ensemble basis:

        |e_a^(m),0>
        |g,1>
        '''
        n_cart_tot = self.n_nuclei * 3
        n_cart_mol = self._n_nuc_mol * 3

        basis_NACs = np.zeros((self._n_dim, self._n_dim, n_cart_tot))

        for i, state_i in enumerate(self._state_pairs):
            mi, ai, ni = state_i

            for j, state_j in enumerate(self._state_pairs):
                mj, aj, nj = state_j

                if mi == -1 or mj == -1:
                    continue

                if mi != mj:
                    continue

                i0 = mi * n_cart_mol
                i1 = (mi + 1) * n_cart_mol

                basis_NACs[i, j, i0:i1] = mol_basis_NACs[mi, ai, aj, :]

        return basis_NACs



