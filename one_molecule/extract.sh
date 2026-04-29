#!/bin/bash
# extract data
# remove old files
rm ex.dat
rm u.dat
rm grad_*.dat
rm nac_*.dat
rm tdp_*.dat
rm gs.dat
# grab energies
grep 'FINAL ENERGY:' output.out | awk '{print $3}' > gs.dat
grep 'Osc. (a.u.)' -A 4 output.out | tail -3 | awk '{print $2}' > ex.dat
# grab dipole moments
grep 'Transition dipole moments:' -A 6 output.out | tail -3  | awk '{print $2 , $3 ,  $4}' > u.dat
# grab gradients
grep 'Total Ground State Gradient (Hartree/Bohr)' -A 13 output.out | tail -10 | awk '{print $2 , $3 , $4}' > grad_0.dat
grep 'Root 1: Total Excited State Gradient (Hartree/Bohr)' output.out -A 13 | tail -10 | awk '{print $2 , $3, $4}' >  grad_1.dat
grep 'Root 2: Total Excited State Gradient (Hartree/Bohr)' output.out -A 13 | tail -10 | awk '{print $2 , $3, $4}' >  grad_2.dat
grep 'Root 3: Total Excited State Gradient (Hartree/Bohr)' output.out -A 13 | tail -10 | awk '{print $2 , $3, $4}' >  grad_3.dat
# grab nacs
grep 'Nonadiabatic coupling between ground state and excited state 1' output.out -A 14 | tail -10 | awk '{print $2 , $3 , $4}' > nac_01.dat
grep 'Nonadiabatic coupling between ground state and excited state 2' output.out -A 14 | tail -10 | awk '{print $2 , $3 , $4}' > nac_02.dat
grep 'Nonadiabatic coupling between ground state and excited state 3' output.out -A 14 | tail -10 | awk '{print $2 , $3 , $4}'> nac_03.dat
grep 'Nonadiabatic coupling between excited state 1 and excited state 2' output.out -A 14 | tail -10 | awk '{print $2 , $3 , $4}' > nac_12.dat
grep 'Nonadiabatic coupling between excited state 1 and excited state 3' output.out -A 14 | tail -10 | awk '{print $2 , $3 , $4}' > nac_13.dat
grep 'Nonadiabatic coupling between excited state 2 and excited state 3' output.out -A 14 | tail -10 | awk '{print $2 , $3 , $4}' > nac_23.dat
# grab transition dipole derivatives
grep '   0 ->    1' output.out -A 40 | grep 'Atom        dMuX/dX' -A 12 | tail -11 | awk '{print $2 , $3 , $4}' > tdp_01_x.dat
grep '   0 ->    1' output.out -A 40 | grep 'Atom        dMuY/dX' -A 12 | tail -11 | awk '{print $2 , $3 , $4}' > tdp_01_y.dat
grep '   0 ->    1' output.out -A 40 | grep 'Atom        dMuZ/dX' -A 12 | tail -11 | awk '{print $2 , $3 , $4}' > tdp_01_z.dat
####################################
grep '   0 ->    2' output.out -A 40 | grep 'Atom        dMuX/dX' -A 12 | tail -11 | awk '{print $2 , $3 , $4}' > tdp_02_x.dat
grep '   0 ->    2' output.out -A 40 | grep 'Atom        dMuY/dX' -A 12 | tail -11 | awk '{print $2 , $3 , $4}' > tdp_02_y.dat
grep '   0 ->    2' output.out -A 40 | grep 'Atom        dMuZ/dX' -A 12 | tail -11 | awk '{print $2 , $3 , $4}' > tdp_02_z.dat
########################################
grep '   0 ->    3' output.out -A 40 | grep 'Atom        dMuX/dX' -A 12 | tail -11 | awk '{print $2 , $3 , $4}' > tdp_03_x.dat
grep '   0 ->    3' output.out -A 40 | grep 'Atom        dMuY/dX' -A 12 | tail -11 | awk '{print $2 , $3 , $4}' > tdp_03_y.dat
grep '   0 ->    3' output.out -A 40 | grep 'Atom        dMuZ/dX' -A 12 | tail -11 | awk '{print $2 , $3 , $4}' > tdp_03_z.dat
####################################
# done

echo 'Done'
