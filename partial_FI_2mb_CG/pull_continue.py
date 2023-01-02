import os
import sys

pull_grompp='gmx grompp -f pull.mdp -o pull_continue.tpr -c pull.gro -p sys_water.top -n index.ndx -maxwarn 1'
pull_run = 'mdrun -deffnm pull_continue -v -nb gpu'
os.system(pull_grompp)
os.system(pull_run)

#dyn_grompp = 'gmx grompp -f dynamic_tension.mdp -o sys_hydro_run.tpr -c pull.gro -p sys_water.top -maxwarn 1'
#dyn_run = 'mdrun -deffnm sys_hydro_run -v -nb gpu'
#os.system(dyn_grompp)
#os.system(dyn_run)
