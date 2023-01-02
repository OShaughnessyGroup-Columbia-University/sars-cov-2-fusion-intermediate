import os, sys
import mdtraj as md


def write_topology(DPPC_num, W_num):

    '''
    proteins = zip(['PROA_P', 'PROB_P', 'PROC_P', 'PROD_P', 'PROE_P', 'PROF_P'], [1,1,1,1,1,1])
    '''
    #proteins = zip(['Protein_A', 'Protein_D'], [3, 3])

    #NA_num = len(table[table['name']=='NA+'])
    #CL_num = len(table[table['name']=='CL-'])
    #molecules = zip(['DPPC','Protein_A', 'Protein_D', 'DPPC', 'W','NA+','CL-'],[DPPC_num1, 3, 3, DPPC_num2, W_num, NA_num, CL_num])
    molecules = zip(['Protein_A', 'DPPC', 'W'],[1, DPPC_num, W_num])

    top = open('sys_dry.top', "w")
    sys.stdout = top
    print ('#include "martini_v2.2.itp"')
    print ('#include "martini_v2.0_solvents.itp"')
    print ('#include "martini_v2.0_ions.itp"')
    print ('#include "martini_v2.0_lipids_all_201506.itp"')

    print ( '#include "Protein_A.itp"')
    print ( '#include "Protein_D.itp"\n')

    #print('; Include Position restraint file\n#ifdef POSRES\n#include "protein_unchanged.itp"\n#endif\n')

    
    title = 'MARTINI SYSTEM'
    print ('[ system ]\n; name\n%s\n\n[ molecules ]\n; name  number'%title)
    print ("\n".join("%-10s %7d"%i for i in molecules))
    top.close()


if __name__ == '__main__':

                
    os.system('source /usr/local/gromacs/bin/GMXRC')
    
    a = md.load('ini2.gro').topology
    table, bonds = a.to_dataframe()
    DPPC_num = len(table[(table['resName']=='DPPC') & (table['name']=='NC3')])
    W_num = len(table[table['resName']=='W'])
    write_topology(DPPC_num, W_num)
    #os.system('gmx grompp -f minimization.mdp -c sys_water.gro -p sys_hydro.top -o ions.tpr -maxwarn 1')
    # select group 3 "water" to fill the ions
    #os.system('echo 3 | gmx genion -s ions.tpr -o sys_hydro.gro -p sys_hydro.top -pname NA+ -nname CL- -neutral -conc 0.15')



