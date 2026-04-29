#!/bin/bash


#    echo "Timestep $j"

for i in {1..1}
    do
        echo "Molecule $i"

        cd traj_$i

	cp step/gs.dat .
        cp step/ex.dat .
        cp step/u.dat .
        cp step/grad_*.dat .
        cp step/nac_*.dat .
        cp step/tdp_*.dat .
        cd ..
done
python data_timestep.py
pwd
echo "DONE"

#    mv energy.npy energy_$j.npy
#    mv gradient.npy gradient_$j.npy
#    mv nac.npy nac_$j.npy
#    mv mu.npy mu_$j.npy

