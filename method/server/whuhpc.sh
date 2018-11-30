#!/bin/bash
module load matlab/R2017b
matlab -nodesktop -nosplash -r "test_parallel_v2('x270', 1)" -logfile test_parallel.log