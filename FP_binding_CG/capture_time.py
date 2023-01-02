import mdtraj as md
import numpy as np
import os

n_runs = 10

capture_time = np.zeros(n_runs)

for i in range (10):
    
    path = 'run%d/'%(i+1)
    os.system('echo "1|13 \n q" | gmx make_ndx -f %ssys_hydro_run.gro -o %soutput.ndx'%(path, path))
    os.system('echo 16 | gmx trjconv -f %ssys_hydro_run.xtc -o %ssys_hydro_run_pbc_promb_10.xtc -skip \
               10 -pbc mol -ur compact -s %ssys_hydro_run.tpr -n %soutput.ndx' %(path, path, path, path))
    os.system('echo 16 | gmx trjconv -f %ssys_hydro_run.gro -o %ssys_hydro_run_promb.gro -n %soutput.ndx \
               -s %ssys_hydro_run.tpr'%(path, path, path, path))

    traj = md.load_xtc('%ssys_hydro_run_pbc_promb_10.xtc'%path, top = '%ssys_hydro_run_promb.gro'%path)
    xyz = traj.xyz
    top = traj.topology
    table, bonds = top.to_dataframe()
    z_dppc = xyz[:,table.name == 'NC3',2]
    
    z_midplane = np.mean(z_dppc, axis = 1)
    z_mb_upper = np.zeros(np.shape(z_midplane))
    z_mb_lower = np.zeros(np.shape(z_midplane))
    z_pro_upper = np.zeros(np.shape(z_midplane))
    z_pro_lower = np.zeros(np.shape(z_midplane))

    for j in range(len(z_midplane)):
        temp = z_dppc[j,:]
        z_mb_upper[j] = np.mean(temp[temp>z_midplane[j]])
        z_mb_lower[j] = np.mean(temp[temp<z_midplane[j]]) 

        temp = xyz[j,table.resName != 'DPPC',2]
        z_pro_upper[j] = np.max(temp)
        z_pro_lower[j] = np.min(temp)

    dt = 1e-3
    t = np.arange(len(z_midplane))*dt
    capture_time[i] = np.max(t[~((z_pro_upper > z_mb_lower) & (z_pro_lower < z_mb_upper))]) + dt

print(capture_time)
