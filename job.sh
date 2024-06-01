#! /bin/bash

#HQ --nodes=1
#HQ --time-request=70m
#HQ --time-limit=70m

##HQ --stdout testdir/stdout.txt
##HQ --stderr testdir/stderr.txt

# Launch model server, send back server URL
# and wait to ensure that HQ won't schedule any more jobs to this allocation.

function get_avaliable_port {
    # Define the range of ports to select from
    MIN_PORT=50000
    MAX_PORT=65535

    # Generate a random port number
    port=$(shuf -i $MIN_PORT-$MAX_PORT -n 1)

    # Check if the port is in use
    while lsof -Pi :$port -sTCP:LISTEN -t >/dev/null; do
        # If the port is in use, generate a new random port number
        port=$(shuf -i $MIN_PORT-$MAX_PORT -n 1)
    done

    echo $port
}

port=$(get_avaliable_port)
export PORT=$port
module load compiler/intel/2022.2 
module load mpi/openmpi/4.1

# Assume that server sets the port according to the environment variable 'PORT'.

#cd /home/hd/hd_hd/hd_uv175/gpfs/hd_uv175-uq_proj/umbridge/models/testmodel ./minimal-server & # CHANGE ME! <---------------------------------------- 24.02.23
python /home/hd/hd_hd/hd_uv175/gpfs/hd_uv175-uq_proj/server_multilevel.py & 

load_balancer_dir="${HOME}/gpfs/hd_uv175-uq_proj/umbridge/hpc" # CHANGE ME!
# echo $HQ_NODE_FILE > "hq_node_file_$HQ_JOB_ID.txt"
# printenv > "env_$HQ_JOB_ID.txt"


host=$(hostname -I | awk '{print $1}')

echo "Waiting for model server to respond at $host:$port..."
while ! curl -s "http://$host:$port/Info" > /dev/null; do
    sleep 1
done
echo "Model server responded"

# Write server URL to file identified by HQ job ID.
mkdir -p "$load_balancer_dir/urls"
echo "http://$host:$port" > "$load_balancer_dir/urls/url-$HQ_JOB_ID.txt"

sleep infinity # keep the job occupied
