clear
addpath(genpath(pwd))
debug_flag = 0;
run_mode = 'parallel';%'parallel'; %'serial';
dataset_name_set = {'CUHK01', 'GRID', 'PRID450s', 'VIPeR', 'CUHK03detected', 'CUHK03labeled'};

for i = 1:length(dataset_name_set)
    test_baseline(debug_flag,run_mode, dataset_name_set{i});
end