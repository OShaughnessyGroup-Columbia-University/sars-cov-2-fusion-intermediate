import os, sys
import mdtraj as md

def read_W_num(file):

    a = md.load('%s'%file).topology
    table, bonds = a.to_dataframe()
    W_num = len(table[table['resName']=='W'])

    return W_num



if __name__ == '__main__':


    for i in range(10):
        
        path = 'run%d/'%(i+1)
        
        os.system('cp martini_files/* %s'%path)
        os.system('cp ../4feb21_ei_viral_mb_final/Protein_*.itp %s'%path)
        
        os.system('gmx grompp -f %sminimization-vaccum.mdp -c %sei_2mb_ini%d.gro -p %ssys_vac.top -o %svac-mini.tpr -maxwarn 1'%(path, path, i+1, path, path))
        os.system('gmx mdrun -deffnm %svac-mini -v'%(path))

        os.system('gmx solvate -cp %svac-mini.gro -cs %swater_martini.gro -o %ssys_water.gro -radius 0.21'%(path, path, path))
        
        W_num = read_W_num('%ssys_water.gro'%path)
        os.system('cp %ssys_vac.top %ssys_water.top'%(path, path))
        os.system('echo "W       %d" >>  %ssys_water.top'%(W_num, path))
        





