#!/bin/bash
cd ../dafe/
module load matlab/R2017b
debug_flag=0
run_mode="'parallel'"
matlab -nodesktop -nosplash -nodisplay -r \
"test_dafe($debug_flag,$run_mode); quit"