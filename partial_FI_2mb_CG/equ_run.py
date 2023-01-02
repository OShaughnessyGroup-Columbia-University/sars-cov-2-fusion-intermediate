import os
import sys

pull_grompp='gmx grompp -f dynamic_tension.mdp -o sys_hydro_run.tpr -c hold.gro -p sys_water.top -n index_3fpnh.ndx -maxwarn 1'
pull_run = 'mdrun -deffnm sys_hydro_run -v -nb gpu'
os.system(pull_grompp)
os.system(pull_run)
