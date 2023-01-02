import os
import mdtraj as md
import sys

def write_topology(file, path):
    
    a = md.load('%s'%file).topology
    table, bonds = a.to_dataframe()
    DPPC_num = len(table[(table['resName']=='DPPC') & (table['name']=='NC3')])
    W_num = len(table[table['resName']=='W'])

    molecules = zip(['Protein_A', 'DPPC', 'W'],[1, DPPC_num, W_num])

    top = open('%ssys_hydro.top'%path, "w")
    sys.stdout = top
    print ('#include "martini_v2.2.itp"')
    print ('#include "martini_v2.0_solvents.itp"')
    print ('#include "martini_v2.0_ions.itp"')
    print ('#include "martini_v2.0_lipids_all_201506.itp"')

    print ( '#include "Protein_A.itp"')
    
    #print('; Include Position restraint file\n#ifdef POSRES\n#include "protein_unchanged.itp"\n#endif\n')

    
    title = 'MARTINI SYSTEM'
    print ('[ system ]\n; name\n%s\n\n[ molecules ]\n; name  number'%title)
    print ("\n".join("%-10s %7d"%i for i in molecules))
    top.close()




i = 0    
path = ''

os.system('echo 0 | gmx editconf -f s2eifpc_cg.pdb -o %sbox%d.pdb -box 7 7 10 -center 3.5 3.5 5 -princ -rotate 0 %d 0'%(path, i+1, i*30))
os.system('python2.7 insane2.py -f %sbox%d.pdb -o %ssys_vac%d.gro -p %ssys_vac.top -pbc rectangular -box 7,7,10 -l DPPC -u DPPC -center -dm 4.5'%(path, i+1, path, i+1, path))


os.system('cp itp_files/* %s'%path)

os.system('gmx grompp -f minimization-vaccum.mdp -c %ssys_vac%d.gro -p %ssys_vac.top -o %svac-mini.tpr -maxwarn 1'%(path, i+1, path, path))
os.system('gmx mdrun -deffnm %svac-mini -v'%(path))    

os.system('gmx solvate -cp %svac-mini.gro -cs water_martini.gro -o %ssys_water.gro -radius 0.21'%(path, path))
    
write_topology('%ssys_water.gro'%path, path)
os.system('gmx grompp -f minimization.mdp -c %ssys_water.gro -p %ssys_hydro.top -o %sions.tpr -maxwarn 1'%(path, path, path))
# select group 3 "water" to fill the ions
os.system('echo 14 | gmx genion -s %sions.tpr -o %ssys_hydro.gro -p %ssys_hydro.top -pname NA+ -nname CL- -neutral -conc 0.15'%(path, path, path))


os.system('gmx grompp -f minimization.mdp -c %ssys_hydro.gro -p %ssys_hydro.top -o %ssys_hydro_min.tpr -maxwarn 2'%(path, path, path))
os.system('gmx mdrun -deffnm %ssys_hydro_min -v'%path)
os.system('gmx grompp -f equilibration.mdp -c %ssys_hydro_min.gro -p %ssys_hydro.top -o %ssys_hydro_eq.tpr'%(path, path, path))
os.system('gmx mdrun -deffnm %ssys_hydro_eq -v'%path)
#os.system('gmx grompp -f dynamic_tension.mdp -c %ssys_hydro_eq.gro -p %ssys_hydro.top -o %ssys_hydro_run.tpr -maxwarn 1'%(path, path, path))


for i in range (10):

    path = 'run%d/'%(i+1)

    os.system('mkdir %s'%path)
    os.system('cp sys_hydro_eq.gro %ssys_hydro_eq.gro'%path)
    os.system('cp sys_hydro.top %ssys_hydro.top'%path)
    os.system('cp dynamic_tension.mdp %sdynamic_tension.mdp'%path)
    os.system('cp itp_files/* %s'%path)
    os.system('cp water_martini.gro %swater_martini.gro'%path)



