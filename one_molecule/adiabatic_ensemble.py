# -*- coding: utf-8 -*-
"""
Created on Sun Mar  8 09:37:46 2026

@author: Sourav
"""

import numpy as np

class AdiabaticStates():
    def __init__(self, n_states, n_nuclei) -> None:
        N = n_states

        self._n_states = N                           # Total size of polaritonic Hamiltonian
        self._n_nuclei = n_nuclei                    # Total number of nuclei
        self._hamiltonian = np.zeros((N, N))         # Hamiltonian    
        self._dH = np.zeros((N, N, n_nuclei*3))      # Hamiltonian derivative
        self._diagonalized = False

        self.eigen_vals = np.zeros(N)               # For storing polaritonic adiabatic state energies
        self.eigen_vecs = np.zeros((N, N))          # For storing polaritonic adiabatic states
        self.NACs = np.zeros((N, N, n_nuclei*3))    # For storing NACs between polaritonic adiabatic states
        self.eigen_val_gradients = np.zeros((N, n_nuclei*3))        # Gradients 
        self.eigen_vec_gradients = np.zeros((N, N, n_nuclei*3))     # Gradients

    @property
    def n_states(self):
        return self._n_states      # returns n_states

    @property
    def n_nuclei(self):
        return self._n_nuclei     # returns n_nulei 

    @property
    def hamiltonian(self):
        return self._hamiltonian   # returns hamiltonian 

    @hamiltonian.setter
    def hamiltonian(self, matrix: np.ndarray):
        if matrix.shape != self._hamiltonian.shape:
            raise ValueError(
                f'Atempting to set Hamiltonian to a size of {matrix.shape} '
                f'when it should be {self._hamiltonian.shape}'
            )
        self._hamiltonian = matrix
        self._hamiltonian.flags['WRITEABLE'] = False

    @property
    def dH(self):
        return self._dH      # returns hamiltonian derivative

    @dH.setter
    def dH(self, matrix: np.ndarray):
        if matrix.shape != self._dH.shape:
            raise ValueError(
                f'Atempting to set Hamiltonian gradient to a size of {matrix.shape} '
                f'when it should be {self._dH.shape}'
            )
        self._dH = matrix
        self._dH.flags['WRITEABLE'] = False

    
    # Initialize with zeros    
    def zero(self):
        self._hamiltonian = np.zeros_like(self._hamiltonian)
        self._dH = np.zeros_like(self._dH)
        self._diagonalized = False

        self.eigen_vals = np.zeros_like(self.eigen_vals)
        self.eigen_vecs = np.zeros_like(self.eigen_vecs)
        self.NACs = np.zeros_like(self.NACs)
        self.eigen_val_gradients = np.zeros_like(self.eigen_val_gradients)
        self.eigen_vec_gradients = np.zeros_like(self.eigen_vec_gradients)

    
    # Compute adiabatic polaritonic state energies and vectors
    def diagonalize_H(self, ref_eig_vecs=None, swap_signs=False):
        e_vals, e_vecs = np.linalg.eigh(self._hamiltonian)    # Diagonalize symmetric matrix
        order = np.argsort(e_vals)                            # Arrange in order of increasing energies

        self.eigen_vals = e_vals[order]
        self.eigen_vecs = e_vecs[:, order]
        self.eigen_vecs[:, 0] *= -1  
        
        if swap_signs:
            for i in range(self.eigen_vecs.shape[1]):
                idx = np.argmax(np.abs(self.eigen_vecs[:, i]))              # check sign of largest value of each column vector
                if self.eigen_vecs[idx, i] < 0:
                    self.eigen_vecs[:, i] *= -1

        if ref_eig_vecs is not None:                                        # check sign with previous step
            corrected_evecs = np.zeros_like(self.eigen_vecs)
            for i in range(self.eigen_vecs.shape[1]):
                dot = np.dot(self.eigen_vecs[:, i], ref_eig_vecs[:, i])
                sgn = 1.0 if dot >= 0 else -1.0
                corrected_evecs[:, i] = self.eigen_vecs[:, i] * sgn
            self.eigen_vecs = corrected_evecs

        self._diagonalized = True
        return self.eigen_vals, self.eigen_vecs                           # return adiabatic states and energies

    
    # Compute adiabatic gradient
    def calc_eigen_value_gradient(self):
        C = self.eigen_vecs
        n_cart = self.dH.shape[2]
        grad = np.zeros((self._n_states, n_cart))
        for k in range(n_cart):
            dHk = self.dH[:, :, k]
            grad[:, k] = np.diag(np.transpose(C) @ dHk @ C)

        self.eigen_val_gradients = grad
        self.eigen_val_gradients.flags['WRITEABLE'] = False
        return self.eigen_val_gradients                                 # return gradients  
 



      # Compute NACs between polaritonic adiabtic states
    def calc_NA_coupling(self, basis_NACs):
          couplings = np.zeros((self._n_states, self._n_states, self._n_nuclei*3))
          for i, E_i in enumerate(self.eigen_vals):
              C_i = self.eigen_vecs[:, i]
              for j, E_j in enumerate(self.eigen_vals):
                  if i == j:
                      continue
                  C_j = self.eigen_vecs[:, j]
                  inverse_energies = 1.0 / (E_j - E_i)
                  term1 = inverse_energies * np.einsum('m,mnk,n->k', C_i, self.dH, C_j)
                  term2 = np.einsum('m,mnk,n->k', C_i, basis_NACs, C_j)
                  couplings[i, j] = term1 + term2

          self.NACs = couplings
          return self.NACs                                               # compute nacs













