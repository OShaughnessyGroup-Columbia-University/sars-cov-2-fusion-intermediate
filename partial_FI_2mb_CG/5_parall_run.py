import os,sys

index = int(sys.argv[1])
path = 'run%d/'%index

os.system('mkdir %s'%path)
os.system('gmx grompp -f dynamic_tension.mdp -o %ssys_hydro_run_continue.tpr -p sys_water.top -c sys_hydro_run.gro -maxwarn 1'%path)
os.system('mdrun -deffnm %ssys_hydro_run_continue -v -nb gpu'%path)
