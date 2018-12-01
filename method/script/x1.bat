cd ..\dafe\
set debug_flag=1
set run_mode="'serial'"
start matlab.exe -nodesktop -nosplash -r ^
"test_dafe(%debug_flag%,%run_mode%); quit"