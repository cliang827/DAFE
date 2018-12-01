cd ..\dafe\
set debug_flag=0
set run_mode="'parallel'"
start matlab.exe -nodesktop -nosplash -r ^
"test_dafe(%debug_flag%,%run_mode%); quit"