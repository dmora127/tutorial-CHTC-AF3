#!/bin/bash

# Run the container with inference script in the background
apptainer exec \
    --scratch /tmp \
    --scratch /var/tmp \
    --workdir "$(pwd)" \
    --pwd "$(pwd)" \
    --bind "$(pwd)" \
    --no-home \
    --nv \
    --containall \
    --env XLA_PYTHON_CLIENT_PREALLOCATE=false \
    alphafold3.minimal.22Jan2025.sif \
    /bin/bash inference_pipeline.sh "$@" \
    1> job.out 2> job.err &
pid=$!

echo "Started inference_pipeline.sh inside container with PID: $pid"

# Monitor GPU usage every 30s until the process exits
while kill -0 "$pid" 2>/dev/null; do
    #echo "===== $(date) ====="
    nvidia-smi --query-compute-apps=pid,process_name,used_memory,timestamp \
           --format=csv,noheader,nounits \
           | tee -a nvidia_raw.log \
           | awk -F, -v pid="$pid" '$1 == pid {print "PID=" $1 " NAME=" $2 " GPU-MEM(MB)=" $3 " Timestamp=" $4}' >> gpu_usage.log
    sleep 10
done

wait $pid

echo "Process $pid finished."
