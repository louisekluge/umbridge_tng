#! /bin/bash

# Note: For runs on systems without SLURM, replace the slurm allocator by
# hq worker start &


hq alloc add slurm --time-limit 630m \
                   --idle-timeout 5m \
                   --backlog 1 \
                   --workers-per-alloc 1 \
                   --max-worker-count 4 \
                   -- -p "single" \
                    --tasks-per-node=64 \
                    --mem=235gb