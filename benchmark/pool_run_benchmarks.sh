#!/bin/bash

# change to directory of the script
cd "$(dirname "$0")"

username="wuerfelh"
all_nodes=$(seq -f pool%g 1 46)
used_nodes=$(ssh ${username}@pool.physik.hu-berlin.de "squeue | awk '/pool/ {print \$8}'")
free_nodes=$(echo "$all_nodes" | grep -Fxv -f <(echo "$used_nodes"))
echo "According to squeue, the follwing nodes are free"
echo $free_nodes

while true; do
    read -p "Where to Benchmark? (type 'check name' to run deep check) " machine

    # Check if the input starts with 'check'
    if [[ $machine == check* ]]; then
        # Extract the name after 'check'
        name=$(echo $machine | awk '{print $2}')

        # Call the checkusage script with the extracted name
        ./pool_checkusage.sh $name
    else
        break
    fi
done

# goto ND root
cd ".."

# set workdir
# create direcoty with current timestamp
dir="/data/scratch/wuerfelh/benchmark_$(date +%Y-%m-%d_%H%M%S)/NetworkDynamics"

echo "Connect to $machine as $username and create $dir"
# create workdid if needed
ssh ${username}@${machine}.physik.hu-berlin.de mkdir -p ${dir}

echo "Copy files to pool: $dir"
rsync -az --recursive --exclude="*.data" \
                      --exclude "*.txt" \
                      --exclude "*.md" \
                      --exclude "*.pdf" --delete . ${username}@${machine}.physik.hu-berlin.de:${dir}

ssh ${username}@${machine}.physik.hu-berlin.de -t "export JULIA_NUM_THREADS=4; cd ${dir}/benchmark; julia --startup-file=no run_benchmarks.jl --no-data-export $@"

# echo "Fetch resulting txt and pdf files"
rsync -a ${username}@${machine}.physik.hu-berlin.de:${dir}/benchmark/*.txt benchmark/
rsync -a ${username}@${machine}.physik.hu-berlin.de:${dir}/benchmark/*.pdf benchmark/
