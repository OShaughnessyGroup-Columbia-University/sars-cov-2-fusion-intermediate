import os, sys
import mdtraj as md
import numpy as np

def read_DPPC_num(path):

    a = md.load('%stop.pdb'%path).topology
    table, bonds = a.to_dataframe()
    DPPC_num1 = len(table[(table['resName']=='DPPC') & (table['name']=='NC3')])

    b = md.load('%sbottom.pdb'%path).topology
    table, bonds = b.to_dataframe()
    DPPC_num2 = len(table[(table['resName']=='DPPC') & (table['name']=='NC3')])

    return DPPC_num1, DPPC_num2

def write_topology(path, DPPC_num1, DPPC_num2):


    molecules = zip(['DPPC','Protein_A', 'Protein_D', 'DPPC'],[DPPC_num1, 3, 3, DPPC_num2])

    top = open('%ssys_vac.top'%path,"w")
    sys.stdout = top
    print ('#include "martini_v2.2.itp"')
    print ('#include "martini_v2.0_solvents.itp"')
    print ('#include "martini_v2.0_ions.itp"')
    print ('#include "martini_v2.0_lipids_all_201506.itp"')

    print ( '#include "Protein_A.itp"')
    print ( '#include "Protein_D.itp"\n')

    
    title = 'MARTINI SYSTEM'
    print ('[ system ]\n; name\n%s\n\n[ molecules ]\n; name  number'%title)
    print ("\n".join("%-10s %7d"%i for i in molecules))
    top.close()


for i in range(10):
    
    path = 'run%d/'%(i+1)
    
    
    os.system('mkdir %s'%path)
    os.system('cp ini%d.gro %s'%(i+1, path))
    traj1 = md.load('%sini%d.gro'%(path, i+1))
    top1 = traj1.topology
    xyz1 = traj1.xyz
    table1, bonds1 = top1.to_dataframe()
    table1.loc[ table1.chainID==0, 'resSeq'] = table1['resSeq'][table1.chainID==0].values + 815
    table1.loc[ table1.chainID==1, 'resSeq'] = table1['resSeq'][table1.chainID==1].values + 815
    table1.loc[ table1.chainID==2, 'resSeq'] = table1['resSeq'][table1.chainID==2].values + 815

    table1.loc[ table1.chainID==3, 'resSeq'] = table1['resSeq'][table1.chainID==3].values + 685
    table1.loc[ table1.chainID==4, 'resSeq'] = table1['resSeq'][table1.chainID==4].values + 685
    table1.loc[ table1.chainID==5, 'resSeq'] = table1['resSeq'][table1.chainID==5].values + 685
    
    table_pro =  table1[table1.chainID < 3]
    xyz_pro = xyz1[:,table1.chainID < 3,:]
    xyz_tm1 = np.mean(xyz_pro[:,(table_pro.resSeq==1214)&(table_pro.name=='BB'),:],axis=1)
    z_tm1 = xyz_tm1[0][2]

    os.system('cp top_mb_hydro.gro %s'%path)

    os.system('echo "q\n" | gmx make_ndx -f %stop_mb_hydro.gro -o %smb_index.ndx'%(path, path))
    os.system('echo 2 | gmx editconf -f %stop_mb_hydro.gro -o %stop_mb.gro -n %smb_index.ndx'%(path, path, path))
    traj2 = md.load('%stop_mb.gro'%path)
    top2 = traj2.topology
    xyz2 = traj2.xyz
    table2, bonds2 = top2.to_dataframe()
    z_dppc2 = np.mean(xyz2[0,table2.resName=='DPPC',2])

    sep = 20+2.5 - (z_dppc2 - z_tm1)
    os.system('gmx editconf -f %stop_mb.gro -o %stop_mb_2.gro -translate 0 0 %f'%(path, path, sep))
    os.system('gmx editconf -f %stop_mb_2.gro -o %stop_temp.pdb'%(path, path))
    os.system('gmx editconf -f %sini%d.gro -o %sbottom_temp.pdb'%(path, i+1, path))
    os.system('tac %stop_temp.pdb | sed "1,2d" | tac > %stop.pdb'%(path, path))
    os.system('tail -n +5 %sbottom_temp.pdb > %sbottom.pdb'%(path, path))
    os.system('cat %stop.pdb %sbottom.pdb > %sboth.pdb'%(path, path, path))
    os.system('gmx editconf -f %sboth.pdb -o %sboth.gro'%(path, path))
    os.system('gmx editconf -f %sboth.gro -bt triclinic -box 31 31 35 -angles 90 90 90 -o %sei_2mb_ini%d.gro'%(path, path, i+1))
    

    DPPC_num1, DPPC_num2 = read_DPPC_num(path)
    write_topology(path, DPPC_num1, DPPC_num2)
