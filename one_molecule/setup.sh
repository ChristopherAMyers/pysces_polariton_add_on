#cp nebpath.xyz path.xyz
#rm temp.xyz
for i in {1..1};     # number of molecules
do
        if [ ! -d "traj_$i" ]; then
        mkdir "traj_$i"
        cp traj_$i.xyz traj_$i
	fi
        cd traj_$i
	cp ../q1.xyz .
################################################################
        echo "Molcule $i"        
        #j=$((k*50*12))      # 10fs = 50 steps (0.2 fs/step)
       # echo $j
       # if [ $j -eq 0 ]; then
       #     head -12 traj_$i.xyz > $k.xyz
       # else
       #     head -$j traj_$i.xyz > temp.xyz
       #     tail -12 temp.xyz > $k.xyz
       # fi
        cat q1.xyz >> traj_$i.xyz
        tail -12 traj_$i.xyz > step.xyz
        cat step.xyz
        mkdir step
        cd step
        cp ../step.xyz .
        cp ../../input.in .
        cd ..
        cd ..
done
