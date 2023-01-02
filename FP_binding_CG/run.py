import os, sys

#os.system('source /usr/local/gromacs/bin/GMXRC')

idx = int(sys.argv[1])
    
path = 'run%d/'%(idx)
os.system('gmx grompp -f %sdynamic_tension.mdp -c %ssys_hydro_eq.gro -p %ssys_hydro.top -o %ssys_hydro_run.tpr -maxwarn 1'%(path, path, path, path))
os.system('mdrun -deffnm %ssys_hydro_run -nb gpu -v'%path) 