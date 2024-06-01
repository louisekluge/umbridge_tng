import numpy as np
import os 
import sys
import h5py
import subprocess
#import mpi4py
import umbridge
import pickle

def write_param_file(WindEn, WindVel, SeedBH, ID):

    doc = os.getcwd()+'param.txt'
    f = open(doc, 'r')
    name = 'param_mod_'+str(ID)+'.txt'
    command = 'touch '+ name
    os.system(command)
    new_path = os.getcwd()+name
    g = open(new_path, 'w')

    for line in f:
        line_s = line.strip()
        columns = line_s.split()
        
        if len(columns)>0:
            if columns[0] == 'OutputDir':
                new_dir = os.getcwd()+"/output_run_" + str(ID) + "/"
                new_string = str(columns[0])+ "\t" + "\t" +"     " + new_dir +"\n"
                line = new_string

        if len(columns)>0:
            if columns[0]=='WindEnergyIn1e51erg':
                new_string = str(columns[0])+ "\t" + "     " + str(WindEn)+"\n"
                line =new_string

        if len(columns)>0:
            if columns[0]=='SeedBlackHoleMass':
                new_string = str(columns[0])+ "\t"  + "\t" + "  "+str(SeedBH)+"\n"
                line =new_string

        if len(columns)>0:
            if columns[0]=='VariableWindVelFactor':
                new_string = str(columns[0])+ "\t" + "     " + str(WindVel)+"\n"
                line =new_string

        g.write(line)

    f.close()
    g.close()
    return(0)

def extract_from_snap(path): # extracts omega_hi for data, omega_star for qoi
    current_path = os.getcwd()
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
    f.close()
    return(omega_hi, omega_star)

def conv_qoi(mstar, volume_mpc):
    # sidelength in Mpc
    HubbleParam = 0.6774
    return(mstar/(27.754*HubbleParam**2*volume_mpc))

def conv_vol20(mstar, z): #Mstar in 10^10 M_sun
    box_z_2 = 9.854221122799048 #in Mpc %a third 
    box_z_1 = 14.77078420600895 #calculated with f['Header'].attrs['Time']*f['Header'].attrs['BoxSize']/0.6774 %half
    box_z_0 = 29.524653085326237
    rho_crit = 1.274 # 10^11 M_sun Mpc^-3
    # sidelength in Mpc
    if z == 0:
        norm = (10*rho_crit*box_z_0**3)
        a = 0.9999999999999996
        return(mstar/norm*a**3)
    if z == 1:
        norm = (10*rho_crit*box_z_1**3)
        a = 0.5002864610575232
        return(mstar/norm*a**3)
    if z == 2:
        norm = (10*rho_crit*box_z_2**3)
        a = 0.33376246942920373
        return(mstar/norm*a**3)

class TNGtest(umbridge.Model):
    def __init__(self):
        super().__init__("forward_test")

    def get_input_sizes(self, config):
        return [3]

    def get_output_sizes(self, config):
        return [1]

    def __call__(self, parameters, config):

        print("Reached call of TNGtest model")

        # Sleep for number of milliseconds defined in env var
        #time.sleep(int(os.getenv("TEST_DELAY", 0)) / 1000)
        run_id = os.environ.get("HQ_JOB_ID")
        print(os.getcwd())
        os.chdir('..')
        print(os.getcwd())
        os.chdir('..')
        print(os.getcwd())
        new_dir = os.getcwd()+'/L25n32/arepo'
        os.chdir(new_dir)
        out_loc = 'mkdir output_run_'+str(run_id)  ### CREATE OUTPUT DIRECTORY
        os.system(out_loc)
        skript = os.getcwd() + "/param_command.py"
        subprocess.run(["python", skript, str(parameters[0][0]), str(parameters[0][1]),str(parameters[0][2]), str(run_id)])
        #os.system(out_loc)
        #subprocess.run(["python", skript, str(parameters[0][0]), str(parameters[0][0]),str(parameters[0][0]), str(run_id)])
        mpi_command = 'mpirun -np 32 ./Arepo param_mod_'+str(run_id)+'.txt &> output_run_'+str(run_id)+'/run_output.txt' ##> output_run_'+str(run_id)+'/run_'+str(run_id)+'.txt'
        os.system(mpi_command)

        current_path = os.getcwd() # do evaluation of outputs
        where = current_path+'/output_run_'+str(run_id)+'/sfr.txt'
        #subprocess.call(["python", "/gpfs/bwfor/work/ws/hd_uv175-uq_proj/L25n32/arepo/param_command.py", parameters[0][0], parameters[0][0], parameters[0][0],run_id])
        #write_param_file(parameters[0][0], parameters[0][0], parameters[0][0], run_id)
        #create = "python3 param_command.py "+str(WindEn)+" "+str(WindVel)+" "+str(SeedBH)+" "+str(run_id)
        #os.system(create)
        #posterior = 2*parameters[0][0]
        test = current_path+'/output_run_'+str(run_id)+'/restartfiles'
        remove_command = 'rm -r '+current_path+'/output_run_'+str(run_id)
        if os.path.exists(test)==False:
            return[[9999999]]
        else:
            data = np.loadtxt(where)    
            sum_totalSm = np.sum(data[:,1])
            return [[sum_totalSm]]

    def supports_evaluate(self):
        return True

class surrogate(umbridge.Model):
    def __init__(self):
        super().__init__("gp_surrogate")

    def get_input_sizes(self, config):
        return [3]

    def get_output_sizes(self, config):
        return [1]

    def __call__(self, parameters, config):

        print("Reached call of surrogate model")

        # Sleep for number of milliseconds defined in env var
        #time.sleep(int(os.getenv("TEST_DELAY", 0)) / 1000)
        run_id = os.environ.get("HQ_JOB_ID")
        print(os.getcwd())
        os.chdir('..')
        print(os.getcwd())
        os.chdir('..')
        print(os.getcwd())
        new_dir = os.getcwd()+'/gp_surrogate'
        os.chdir(new_dir)

        filename = '/home/hd/hd_hd/hd_uv175/gpfs/hd_uv175-uq_proj/gp_surrogate/gp_multioutput.pkl'

        gp = pickle.load(open(filename, 'rb'))

        X_eval = np.array([[parameters[0][0], parameters[0][1], parameters[0][2]]])

        output = gp.predict(X_eval, return_std = True)

        #with open('evaluations.txt', "a") as f:
    
            #f.write(str(parameters[0][0])+'\t' +str(parameters[0][1])+'\t' +str(parameters[0][2])+'\t' + str(output[0]) +'\t' + str(std[0]) + "\n")

        #f.close()
        return[[conv_qoi(output[0], 20.0)]]
        #return [[conv_qoi(output[0], 20.0)]] #only QoI

    def supports_evaluate(self):
        return True

class TNGn32(umbridge.Model):

    def __init__(self):
        self.name = "TNGn32"
        super().__init__("forward_level1")

    def get_input_sizes(self, config):
        return [3]

    def get_output_sizes(self, config):
        return [3]

    def __call__(self, parameters, config):

        print("Reached call of TNGn32 model")

        #hand parameters to script for parameter file preparation
        run_id = os.environ.get("HQ_JOB_ID")
        os.chdir('..')
        os.chdir('..')
        new_dir = os.getcwd()+'/L10n32/arepo'
        os.chdir(new_dir)

        out_loc = 'mkdir output_run_'+str(run_id)  ### CREATE OUTPUT DIRECTORY
        os.system(out_loc)
        skript = os.getcwd() + "/param_command.py"
        subprocess.run(["python", skript, str(np.exp(parameters[0][0])), str(np.exp(parameters[0][1])),str(np.exp(parameters[0][2])), str(run_id)])

        mpi_command = 'mpirun -np 32 ./Arepo param_mod_'+str(run_id)+'.txt &> output_run_'+str(run_id)+'/run_output.txt' ##> output_run_'+str(run_id)+'/run_'+str(run_id)+'.txt'
        os.system(mpi_command)
        current_path = os.getcwd() # do evaluation of outputs
 
        current_path = os.getcwd() # do evaluation of outputs
        where = current_path+'/output_run_'+str(run_id)+'/sfr.txt'
        log_skript = os.getcwd() + "/log_runs_hq.py"
        test = current_path+'/output_run_'+str(run_id)+'/snap_002.hdf5'
        remove_command = 'rm -r '+current_path+'/output_run_'+str(run_id)
        if os.path.exists(test)==False:
            os.system(remove_command)
            return[[9999999, 9999999, 9999999]]
        else:
            hi_z2, star_z2 = extract_from_snap('/output_run_'+str(run_id)+'/snap_000.hdf5')
            hi_z1, star_z1 = extract_from_snap('/output_run_'+str(run_id)+'/snap_001.hdf5')
            hi_z0, star_z0 = extract_from_snap('/output_run_'+str(run_id)+'/snap_002.hdf5')

            subprocess.run(["python", log_skript, str(np.exp(parameters[0][0])), str(np.exp(parameters[0][1])),str(np.exp(parameters[0][2])), str(run_id)])
            os.system(remove_command)
            return [[hi_z0, hi_z1, hi_z2]]
    
    def supports_evaluate(self):
        return True

    def supports_gradient(self):
        return False

class TNGn64(umbridge.Model):

    def __init__(self):
        self.name = "TNGn64"
        super().__init__("forward_level2")

    def get_input_sizes(self, config):
        return [3]

    def get_output_sizes(self, config):
        return [3]

    def __call__(self, parameters, config):
        print("Reached call of TNGn64 model")
        run_id = os.environ.get("HQ_JOB_ID")
        os.chdir('..')
        os.chdir('..')
        new_dir = os.getcwd()+'/L25n64/arepo'
        os.chdir(new_dir)

        out_loc = 'mkdir output_run_'+str(run_id)  ### CREATE OUTPUT DIRECTORY
        os.system(out_loc)
        skript = os.getcwd() + "/param_command.py"
        subprocess.run(["python", skript, str(np.exp(parameters[0][0])), str(np.exp(parameters[0][1])),str(np.exp(parameters[0][2])), str(run_id)])

        mpi_command = 'mpirun -np 64 ./Arepo param_mod_'+str(run_id)+'.txt &> output_run_'+str(run_id)+'/run_output.txt' ##> output_run_'+str(run_id)+'/run_'+str(run_id)+'.txt'

        os.system(mpi_command)
        current_path = os.getcwd() # do evaluation of outputs
        where = current_path+'/output_run_'+str(run_id)+'/sfr.txt'
        log_skript = os.getcwd() + "/log_runs_hq.py"
        test = current_path+'/output_run_'+str(run_id)+'/snap_002.hdf5'
        remove_command = 'rm -r '+current_path+'/output_run_'+str(run_id)
        if os.path.exists(test)==False:
            #os.system(remove_command)
            return[[9999999, 9999999, 9999999]]
        else:
            hi_z2, star_z2 = extract_from_snap('/output_run_'+str(run_id)+'/snap_000.hdf5')
            hi_z1, star_z1 = extract_from_snap('/output_run_'+str(run_id)+'/snap_001.hdf5')
            hi_z0, star_z0 = extract_from_snap('/output_run_'+str(run_id)+'/snap_002.hdf5')

            subprocess.run(["python", log_skript, str(np.exp(parameters[0][0])), str(np.exp(parameters[0][1])),str(np.exp(parameters[0][2])), str(run_id)])
            os.system(remove_command)
            return [[hi_z0, hi_z1, hi_z2]]
    
    def supports_evaluate(self):
        return True

    
class TNGn128(umbridge.Model):

    def __init__(self):
        self.name = "TNGn128"
        super().__init__("forward_level3")

    def get_input_sizes(self, config):
        return [3]

    def get_output_sizes(self, config):
        return [3]

    def __call__(self, parameters, config):
        print("Reached call of TNGn128 model")

        run_id = os.environ.get("HQ_JOB_ID")
        os.chdir('..')
        os.chdir('..')
        new_dir = os.getcwd()+'/L25n128/arepo'
        os.chdir(new_dir)

        out_loc = 'mkdir output_run_'+str(run_id)  ### CREATE OUTPUT DIRECTORY
        os.system(out_loc)
        skript = os.getcwd() + "/param_command.py"
        subprocess.run(["python", skript, str(np.exp(parameters[0][0])), str(np.exp(parameters[0][1])),str(np.exp(parameters[0][2])), str(run_id)])

        mpi_command = 'mpirun -np 256 ./Arepo param_mod_'+str(run_id)+'.txt &> output_run_'+str(run_id)+'/run_output.txt' ##> output_run_'+str(run_id)+'/run_'+str(run_id)+'.txt'
        #--node-list=$HQ_NODE_FILE 
        os.system(mpi_command)

        current_path = os.getcwd() # do evaluation of outputs
        where = current_path+'/output_run_'+str(run_id)+'/sfr.txt'
        log_skript = os.getcwd() + "/log_runs_hq.py"
        test = current_path+'/output_run_'+str(run_id)+'/restartfiles'
        remove_command = 'rm -r '+current_path+'/output_run_'+str(run_id)
        if os.path.exists(test)==False:
            #os.system(remove_command)
            return[[9999999, 9999999, 9999999]]
        else:
            hi_z2, star_z2 = extract_from_snap('/output_run_'+str(run_id)+'/snap_000.hdf5')
            hi_z1, star_z1 = extract_from_snap('/output_run_'+str(run_id)+'/snap_001.hdf5')
            hi_z0, star_z0 = extract_from_snap('/output_run_'+str(run_id)+'/snap_002.hdf5')

            subprocess.run(["python", log_skript, str(np.exp(parameters[0][0])), str(np.exp(parameters[0][1])),str(np.exp(parameters[0][2])), str(run_id)])
            #os.system(remove_command)
            return [[hi_z0, hi_z1, hi_z2]]
    
    def supports_evaluate(self):
        return True
    
class TNGn256(umbridge.Model):

    def __init__(self):
        self.name = "TNGn256"
        super().__init__("forward_level4")

    def get_input_sizes(self, config):
        return [3]

    def get_output_sizes(self, config):
        return [3]

    def __call__(self, parameters, config):
        print("Reached call of TNGn256 model")

        run_id = os.environ.get("HQ_JOB_ID")
        os.chdir('..')
        os.chdir('..')
        new_dir = os.getcwd()+'/L25n256/arepo'
        os.chdir(new_dir)

        out_loc = 'mkdir output_run_'+str(run_id)  ### CREATE OUTPUT DIRECTORY
        os.system(out_loc)
        skript = os.getcwd() + "/param_command.py"
        subprocess.run(["python", skript, str(np.exp(parameters[0][0])), str(np.exp(parameters[0][1])),str(np.exp(parameters[0][2])), str(run_id)])

        mpi_command = 'mpirun -np 512 ./Arepo param_mod_'+str(run_id)+'.txt &> output_run_'+str(run_id)+'/run_output.txt' ##> output_run_'+str(run_id)+'/run_'+str(run_id)+'.txt'
        os.system(mpi_command)

        current_path = os.getcwd() # do evaluation of outputs
        where = current_path+'/output_run_'+str(run_id)+'/sfr.txt'
        log_skript = os.getcwd() + "/log_runs_hq.py"
        test = current_path+'/output_run_'+str(run_id)+'/restartfiles'
        remove_command = 'rm -r '+current_path+'/output_run_'+str(run_id)
        if os.path.exists(test)==False:
            os.system(remove_command)
            return[[9999999]]
        else:
            hi_z2, star_z2 = extract_from_snap('/output_run_'+str(run_id)+'/snap_000.hdf5')
            hi_z1, star_z1 = extract_from_snap('/output_run_'+str(run_id)+'/snap_001.hdf5')
            hi_z0, star_z0 = extract_from_snap('/output_run_'+str(run_id)+'/snap_002.hdf5')

            subprocess.run(["python", log_skript, str(np.exp(parameters[0][0])), str(np.exp(parameters[0][1])),str(np.exp(parameters[0][2])), str(run_id)])
            os.system(remove_command)
            return [[hi_z0, hi_z1, hi_z2]]
    
    def supports_evaluate(self):
        return True

testmodel = TNGtest()
surrogate_gp = surrogate()
model_level1 = TNGn32()
model_level2 = TNGn64()
model_level3 = TNGn128()
model_level4 = TNGn256()


if __name__ == "__main__":
    port = int(os.environ["PORT"])
    print(f"Running TNG server on port {port}.")
    umbridge.serve_models([surrogate_gp, model_level1, model_level2, model_level3, model_level4], port)