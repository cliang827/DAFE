#!/bin/bash
cd ../dafe/
debug_flag=1
run_mode="'serial'"
matlab -nodesktop -nosplash -nodisplay -r \
"test_dafe($debug_flag, $run_mode); quit"