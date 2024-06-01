import numpy as np
import sys 
import os
import h5py

current_path = os.getcwd()

def extract_from_snap(path): # extracts omega_hi for data, omega_star for qoi
    snap_loc = current_path + path
    f = h5py.File(snap_loc)
    boxsize = (f['Header'].attrs['BoxSize']*(10**(-3))/0.6774)**3 # in cMpc^3
    rho_crit = 1.274 # 10^11 M_sun Mpc^-3 AT REDSHIFT z = 0

    masses = np.array(f['PartType0']['Masses'])
    h1_abundance = np.array(f['PartType0']['NeutralHydrogenAbundance'])
    m_hi = np.sum(masses*h1_abundance)*0.75/0.6774/10 # in 10^11 M_sol
    omega_hi = m_hi/boxsize/rho_crit

    m_star = 0
    if len(f.keys())>7:
        m_star = np.sum(f['PartType4']['Masses'])/0.6774/10
    omega_star = m_star/boxsize/rho_crit
    return(omega_hi, omega_star)

par1 = float(sys.argv[1])
par2 = float(sys.argv[2])
par3 = float(sys.argv[3])
run_id = int(sys.argv[4])

where = '/output_run_'+str(run_id)+'/sfr.txt'
here = current_path + where
data = np.loadtxt(here)

sum_totalSm = np.sum(data[:,1])/0.6774
sum_totSfrRate = np.sum(data[:,2])/0.6774
sum_sfrMsunPerYr = np.sum(data[:,3])/0.6774
last_totSMS = data[-1,4]/0.6774
sum_totSMS = np.sum(data[:,4])/0.6774
last_cumMS = data[-1,5]/0.6774

hi_z2, star_z2 = extract_from_snap('/output_run_'+str(run_id)+'/snap_000.hdf5')
hi_z1, star_z1 = extract_from_snap('/output_run_'+str(run_id)+'/snap_001.hdf5')
hi_z0, star_z0 = extract_from_snap('/output_run_'+str(run_id)+'/snap_002.hdf5')

log_path = current_path + "/hq_run64_log.txt"
with open(log_path, "a") as f:

    f.write(str(run_id) + "\t" + str(hi_z2) + "\t"+ str(hi_z1)+ "\t" + str(hi_z0)+ "\t" + str(star_z2)+ "\t"  + str(star_z1) + "\t" + str(star_z0) + "\t" + str(sum_totalSm) + "\t" + str(par1)+ "\t" + str(par2) + "\t" + str(par3)+ "\t" 'n64' +"\n")

f.close()
