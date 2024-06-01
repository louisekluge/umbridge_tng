import os
import sys



def write_param_file(WindEn, WindVel, SeedBH, ID):
    f = open('param.txt', 'r')
    name = 'param_mod_'+str(ID)+'.txt'
    command = 'touch '+ name
    os.system(command)
    g = open(name, 'w')

    for line in f:
        line_s = line.strip()
        columns = line_s.split()
        
        if len(columns)>0:
            if columns[0] == 'OutputDir':
                new_dir = "./output_run_" + str(ID) + "/"
                new_string = str(columns[0])+ "\t" + "     " + new_dir +"\n"
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

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python3 param_command.py <WindEn> <WindVel> <SeedBH> <run_id>")
        sys.exit(1)

    num1 = float(sys.argv[1])
    num2 = float(sys.argv[2])
    num3 = float(sys.argv[3])
    run_id = int(sys.argv[4])

    result = write_param_file(num1, num2, num3, run_id)
    output = "Param file with WindEnergyIn1e51erg "+str(num1)+", VariableWindVelFactor "+str(num2)+", SeedBlackHoleMass "+str(num3)+" generated."
    print(output)
