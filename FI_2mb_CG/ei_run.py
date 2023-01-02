import os,sys

index = int(sys.argv[1])
path = 'run%d/'%index

os.system('gmx grompp -f minimization.mdp -c %ssys_water.gro -p %ssys_water.top -o %sions.tpr -maxwarn 1'%(path, path, path))
# select group 3 "water" to fill the ions
os.system('echo 3 | gmx genion -s %sions.tpr -o %ssys_hydro.gro -p %ssys_water.top -pname NA+ -nname CL- -neutral -conc 0.15'%(path, path, path))
os.system('cp %ssys_water.top %ssys_hydro.top'%(path, path))

os.system('gmx grompp -f minimization.mdp -c %ssys_hydro.gro -p %ssys_hydro.top -o %ssys_hydro_min.tpr -maxwarn 2'%(path, path, path))
os.system('mdrun -deffnm %ssys_hydro_min -v'%path)
os.system('gmx grompp -f equilibration.mdp -c %ssys_hydro_min.gro -p %ssys_hydro.top -o %ssys_hydro_eq.tpr'%(path, path, path))
os.system('mdrun -deffnm %ssys_hydro_eq -v'%path)
os.system('gmx grompp -f dynamic_tension.mdp -c %ssys_hydro_eq.gro -p %ssys_hydro.top -o %ssys_hydro_run.tpr -maxwarn 1'%(path, path, path))
os.system('mdrun -deffnm %ssys_hydro_run -nb gpu -v'%path) 
