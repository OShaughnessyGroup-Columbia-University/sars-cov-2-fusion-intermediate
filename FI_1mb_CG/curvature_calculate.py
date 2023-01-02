import mdtraj as md
import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits import mplot3d
import statsmodels.api as sm
from matplotlib import animation
from IPython.display import HTML
import matplotlib
from scipy.interpolate import interp1d
from scipy import interpolate

trunc = 4.0
s_num = 1001
poly = 4

def read_traj(index):
    
    global traj, xyz_pro, xyz_dppc, xyz_pro_fit_center, table_pro, table
    
    traj = md.load_xtc('run%d/sys_hydro_run_10.xtc'%index, top = 'initial_promb.pdb')
    xyz = traj.xyz
    top = traj.topology
    table, bonds = top.to_dataframe()
    
    table.loc[ table.chainID==0, 'resSeq'] = table['resSeq'][table.chainID==0].values + 815
    table.loc[ table.chainID==1, 'resSeq'] = table['resSeq'][table.chainID==1].values + 815
    table.loc[ table.chainID==2, 'resSeq'] = table['resSeq'][table.chainID==2].values + 815

    table.loc[ table.chainID==3, 'resSeq'] = table['resSeq'][table.chainID==3].values + 685
    table.loc[ table.chainID==4, 'resSeq'] = table['resSeq'][table.chainID==4].values + 685
    table.loc[ table.chainID==5, 'resSeq'] = table['resSeq'][table.chainID==5].values + 685
    
    table_pro =  table[table.chainID < 3]
    xyz_pro = xyz[:,table.chainID < 3,:]
    
    xyz_pro = deal_pbc()
    
    xyz_dppc = xyz[:,(table.resName == 'DPPC') ,:]
    table_pro_fit = table_pro[(table_pro.resSeq > 911)&(table_pro.resSeq < 1238)&(table_pro.name == 'BB')]
    xyz_pro_fit = xyz_pro[:,(table_pro.resSeq > 911)&(table_pro.resSeq < 1238)&(table_pro.name == 'BB'),:]
    chain_length_pro_fit = int(np.shape(xyz_pro_fit)[1]/3)
    xyz_pro_fit1 = xyz_pro_fit[:,:chain_length_pro_fit,:]
    xyz_pro_fit2 = xyz_pro_fit[:,chain_length_pro_fit:2*chain_length_pro_fit,:]
    xyz_pro_fit3 = xyz_pro_fit[:,2*chain_length_pro_fit:,:]
    xyz_pro_fit_center = (xyz_pro_fit1+xyz_pro_fit2+xyz_pro_fit3)/3


def deal_pbc():
    global xyz_pro, table_pro, traj
    
    chain1_center = np.mean(xyz_pro[:,table_pro.chainID==0,:],axis=1)
    chain2_center = np.mean(xyz_pro[:,table_pro.chainID==1,:],axis=1)
    chain3_center = np.mean(xyz_pro[:,table_pro.chainID==2,:],axis=1)

    for j in range(3):
        pbc_cross_1 = np.where(abs(chain1_center[:,j] - chain2_center[:,j])>4)
        for i in range(len(pbc_cross_1[0])):
            if(chain2_center[pbc_cross_1[0][i],j] - chain1_center[pbc_cross_1[0][i],j] > 0):
                xyz_pro[pbc_cross_1[0][i],table_pro.chainID==1,j] = xyz_pro[pbc_cross_1[0][i],table_pro.chainID==1,j] - traj.unitcell_lengths[pbc_cross_1[0][i],j]
            else:
                xyz_pro[pbc_cross_1[0][i],table_pro.chainID==1,j] = xyz_pro[pbc_cross_1[0][i],table_pro.chainID==1,j] + traj.unitcell_lengths[pbc_cross_1[0][i],j]

        pbc_cross_2 = np.where(abs(chain1_center[:,j] - chain3_center[:,j])>4)
        for i in range(len(pbc_cross_2[0])):
            if(chain3_center[pbc_cross_2[0][i],j] - chain1_center[pbc_cross_2[0][i],j] > 0):
                xyz_pro[pbc_cross_2[0][i],table_pro.chainID==2,j] = xyz_pro[pbc_cross_2[0][i],table_pro.chainID==2,j] - traj.unitcell_lengths[pbc_cross_2[0][i],j]
            else:
                xyz_pro[pbc_cross_2[0][i],table_pro.chainID==2,j] = xyz_pro[pbc_cross_2[0][i],table_pro.chainID==2,j] + traj.unitcell_lengths[pbc_cross_2[0][i],j]
    
    return xyz_pro


def fit_curve(frame_num,s_num,trunc,poly):
    
    n_diff = 260
    n_overlap = 20
    
    global xyz_pro_fit_center
    
    x_fit = xyz_pro_fit_center[frame_num,:,0]
    y_fit = xyz_pro_fit_center[frame_num,:,1]
    z_fit = xyz_pro_fit_center[frame_num,:,2]

    x_fit_norm = x_fit - np.mean(x_fit)
    y_fit_norm = y_fit - np.mean(y_fit)
    z_fit_norm = z_fit - np.mean(z_fit)

    x_fit_norm1 = x_fit_norm[:n_diff+n_overlap]
    y_fit_norm1 = y_fit_norm[:n_diff+n_overlap]
    z_fit_norm1 = z_fit_norm[:n_diff+n_overlap]
    points_data1 = np.array([x_fit_norm1,y_fit_norm1,z_fit_norm1])
    covMatrix1 = np.cov(points_data1)
    a1,b1 = np.linalg.eigh(covMatrix1)

    x_rotate1 = np.dot(b1.T, points_data1)[0,:]
    y_rotate1 = np.dot(b1.T, points_data1)[1,:]
    z_rotate1 = np.dot(b1.T, points_data1)[2,:]

    lowess = sm.nonparametric.lowess
    xz1 = lowess(x_rotate1, z_rotate1, frac = 0.1, it=10)
    yz1 = lowess(y_rotate1, z_rotate1, frac = 0.1, it=10)
    xyz1 = [xz1[np.unique(xz1[:,0], return_index=True)[1],1],yz1[np.unique(xz1[:,0], return_index=True)[1],1],xz1[np.unique(xz1[:,0], return_index=True)[1],0]]

    xyz1_final = np.dot(np.linalg.inv(b1.T), xyz1)
    if(np.dot(xyz1_final[:,-1]-xyz1_final[:,0], points_data1[:,-1]-points_data1[:,0])<0):
        xyz1_final = xyz1_final[:, ::-1]


    x_fit_norm2 = x_fit_norm[n_diff-n_overlap:]
    y_fit_norm2 = y_fit_norm[n_diff-n_overlap:]
    z_fit_norm2 = z_fit_norm[n_diff-n_overlap:]
    points_data2 = np.array([x_fit_norm2,y_fit_norm2,z_fit_norm2])

    covMatrix2 = np.cov(points_data2)
    a2,b2 = np.linalg.eigh(covMatrix2)

    x_rotate2 = np.dot(b2.T, points_data2)[0,:]
    y_rotate2 = np.dot(b2.T, points_data2)[1,:]
    z_rotate2 = np.dot(b2.T, points_data2)[2,:]

    lowess = sm.nonparametric.lowess
    xz2 = lowess(x_rotate2, z_rotate2, frac = 0.1, it=10)
    yz2 = lowess(y_rotate2, z_rotate2, frac = 0.1, it=10)
    xyz2 = [xz2[np.unique(xz2[:,0], return_index=True)[1],1],yz2[np.unique(xz2[:,0], return_index=True)[1],1],xz2[np.unique(xz2[:,0], return_index=True)[1],0]]

    xyz2_final = np.dot(np.linalg.inv(b2.T), xyz2)
    if(np.dot(xyz2_final[:,-1]-xyz2_final[:,0], points_data2[:,-1]-points_data2[:,0])<0):
        xyz2_final = xyz2_final[:, ::-1]

    xyz_final = np.append(np.append(xyz1_final[:,:-2*n_overlap], (xyz1_final[:,-2*n_overlap:]+xyz2_final[:,:2*n_overlap])/2 ,axis=1), xyz2_final[:,2*n_overlap:],axis=1)

    tck,u = interpolate.splprep(xyz_final,k=poly,s=trunc)
    u = np.linspace(0,1,s_num)
    fit_out = np.array(interpolate.splev(u,tck))
    fit_out[0,:] = fit_out[0,:] + np.mean(x_fit)
    fit_out[1,:] = fit_out[1,:] + np.mean(y_fit)
    fit_out[2,:] = fit_out[2,:] + np.mean(z_fit)
    
    z_tm_n = xyz_pro_fit_center[frame_num,1213-912,2]
    #out = fit_out[:,fit_out[2] > z_tm_n]
    out = fit_out.copy()
    
    return out[0], out[1], out[2]

def cal_curvature(frame_num,s_num):
    
    global xyz_pro_fit_center
    
    x_inter, y_inter, z_inter = fit_curve(frame_num,s_num,trunc,poly)
    z_tm_n = xyz_pro_fit_center[frame_num,1213-912,2]
    
    s_inter = np.zeros(np.shape(z_inter))
    s_inter[1:] = np.cumsum(np.sqrt((x_inter[1:]-x_inter[:-1])**2 + (y_inter[1:]-y_inter[:-1])**2 + (z_inter[1:]-z_inter[:-1])**2))
    
    s_num_total = round((s_num-1)/(len(z_inter[z_inter > z_tm_n])-1)*(len(z_inter)-1)+1)
    
    s_ei = np.linspace(0, s_inter[-1], s_num_total)
    fsx = interp1d(s_inter, x_inter, kind='cubic')
    fsy = interp1d(s_inter, y_inter, kind='cubic')
    fsz = interp1d(s_inter, z_inter, kind='cubic')
    x_s = fsx(s_ei)
    y_s = fsy(s_ei)
    z_s = fsz(s_ei)
    
    d2rds2 = np.zeros((3, len(s_ei)))
    ds = s_ei[1] - s_ei[0]
    d2rds2[0,1:-1] = (x_s[2:] + x_s[:-2] - 2 * x_s[1:-1])/ds**2
    d2rds2[0,0] = (2*x_s[0] -5*x_s[1] +4*x_s[2] -1*x_s[3])/ds**2
    d2rds2[0,-1] = (2*x_s[-1] -5*x_s[-2] +4*x_s[-3] -1*x_s[-4])/ds**2
    d2rds2[1,1:-1] = (y_s[2:] + y_s[:-2] - 2 * y_s[1:-1])/ds**2
    d2rds2[1,0] = (2*y_s[0] -5*y_s[1] +4*y_s[2] -1*y_s[3])/ds**2
    d2rds2[1,-1] = (2*y_s[-1] -5*y_s[-2] +4*y_s[-3] -1*y_s[-4])/ds**2
    d2rds2[2,1:-1] = (z_s[2:] + z_s[:-2] - 2 * z_s[1:-1])/ds**2
    d2rds2[2,0] = (2*z_s[0] -5*z_s[1] +4*z_s[2] -1*z_s[3])/ds**2
    d2rds2[2,-1] = (2*z_s[-1] -5*z_s[-2] +4*z_s[-3] -1*z_s[-4])/ds**2
    
    curvature_s = np.linalg.norm(d2rds2, axis=0)
    xyz_s = [x_s, y_s, z_s]
    
    return s_ei, curvature_s, xyz_s

def find_residue_position(frame_num):
    
    x_inter, y_inter, z_inter = fit_curve(frame_num,s_num,trunc,poly)
    
    s_inter = np.zeros(np.shape(z_inter))
    s_inter[1:] = np.cumsum(np.sqrt((x_inter[1:]-x_inter[:-1])**2 + (y_inter[1:]-y_inter[:-1])**2 + (z_inter[1:]-z_inter[:-1])**2))
    
    position = np.zeros(1237-911)
    
    for i in range(1237-911):
        position[i] = s_inter[np.argmin(np.linalg.norm([x_inter-xyz_pro_fit_center[frame_num,i,0], y_inter-xyz_pro_fit_center[frame_num,i,1], z_inter-xyz_pro_fit_center[frame_num,i,2]],axis=0))]
    
    return position


if __name__ == '__main__':
    
    #trunc = 1
    #s_num = 200
    #poly = 4


    frame_total = 2000
    index_array = [1,2,3,4,5]
    curvature_t_5 =  np.zeros((len(index_array), frame_total, s_num))
    s_t_5 = np.zeros(np.shape(curvature_t_5))
    resi_s_5 = np.zeros((len(index_array),frame_total,1213-911))

    for j in range(len(index_array)):

        frame_num = 0
        read_traj(index_array[j])

        curvature_t = np.zeros((frame_total, s_num))
        s_t = np.zeros(np.shape(curvature_t))
        #xyz_t = np.zeros((frame_total,3,s_num))

        s_temp, curvature_temp, xyz_temp = cal_curvature(frame_num, s_num)

        curvature_t[0,:] = curvature_temp[:s_num]
        s_t[0,:] = s_temp[:s_num]
        #xyz_t[0,:,:] = xyz_temp[:s_num]

        pos_temp = find_residue_position(frame_num)
        resi_s_5[j,0,:] = pos_temp[:1213-911]

        for i in range (np.shape(curvature_t)[0]-1):
            frame_num = frame_num + 1

            s_temp, curvature_temp, xyz_temp = cal_curvature(frame_num, s_num)
            s_t[i+1,:] = s_temp[:s_num]
            curvature_t[i+1,:] = curvature_temp[:s_num]
            #xyz_t[i+1,:,:] = xyz_temp[:s_num]
            pos_temp = find_residue_position(frame_num)
            resi_s_5[j,i+1,:] = pos_temp[:1213-911]

        curvature_t_5[j,:,:] = curvature_t
        s_t_5[j,:,:] = s_t

    np.save('curvature_trunc_%d.npy'%(int(trunc)),curvature_t_5)
    np.save('arclength_trunc_%d.npy'%(int(trunc)),s_t_5)
    np.save('resipos_trunc_%d.npy'%(int(trunc)),resi_s_5)



