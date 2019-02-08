function test_dafe(debug_flag, run_mode)

debug_flag = 1; 
run_mode = 'serial'; %'parallel'; %'serial';

clearvars -except debug_flag run_mode
clc
close all

init_environment;
init_parameters;


diary([ctrl_para.dir_info.log_file(1:end-4),'.txt']); diary on;

fprintf(1, 'machine_type=%s, run_mode=%s, debug_flag=%d, version=%s\n', ...
    machine_type, run_mode, debug_flag, version_str{1});

fprintf(1, 'dataset=%s, probe_set_num=%d, trial_num=%d, tot_query_times=%d\n', ...
    dataset_name, curr_dataset.probe_set_num, trial_num, ctrl_para.exp.tot_query_times);

fprintf(1, 'include_groundtruth_flag=%d, v_sum_constraint_flag=%d\n\n', ...
    ctrl_para.exp.include_groundtruth_flag, ctrl_para.exp.v_sum_constraint);

clearvars -except run_mode para_test_set dataset_set ctrl_para_set eval_para 

%%
result_file = eval_para.result_file;
para_test_num = size(para_test_set,1);
reid_score = cell(para_test_num,1);
auc_score = cell(para_test_num,1);
difficulty_score = cell(para_test_num,1);
feedback_id = cell(para_test_num,1);
time_result = cell(para_test_num,1);
switch run_mode
    case 'serial'
        for i=1:para_test_num

            fprintf(1,'[%4d/%d] fb_method=''%s'', log(alpha)=%.2f, beta_percentage=%.2f%%, log(gamma)=%.2f, delta=%.2f, fb_num=%d, t=%d\n', ...
                i, para_test_num, ctrl_para_set{i}.model.fb_method, log10(ctrl_para_set{i}.model.alpha), ...
                100*ctrl_para_set{i}.model.beta_percentage, log10(ctrl_para_set{i}.model.gamma), ...
                ctrl_para_set{i}.model.delta, ctrl_para_set{i}.model.fb_num, ctrl_para_set{i}.exp.trial);

            [reid_score{i}, difficulty_score{i}, auc_score{i}, feedback_id{i}, time_result{i}] = ...
                dafe(dataset_set{i}, ctrl_para_set{i}); 
        end
        save(result_file, 'para_test_set', 'eval_para', 'reid_score', 'difficulty_score', ...
            'auc_score', 'feedback_id', 'time_result', '-v7.3');
        [time_in_total, time_each_probe] = analyze_parameters(result_file);
        fprintf(1, '\n time_in_total=%.0f sec., time_each_probe:%.2f sec\n', ...
            time_in_total, time_each_probe);
        
    case 'parallel'
        t_start =tic;
        parfor i=1:para_test_num
            [reid_score{i}, difficulty_score{i}, auc_score{i}, feedback_id{i}, time_result{i}] = ...
                dafe(dataset_set{i}, ctrl_para_set{i});
        end
        if ~isempty(gcp('nocreate'))>0
            delete(gcp('nocreate'))
        end
        time_in_parallel = toc(t_start);
        save(result_file, 'para_test_set', 'eval_para', 'reid_score', 'difficulty_score', ...
            'auc_score', 'feedback_id', 'time_result', 'time_in_parallel', '-v7.3');
        [time_in_serial, time_each_probe] = analyze_parameters(result_file);
        fprintf(1, '\n time_in_parallel=%.0f sec, time_in_total_serial=%.0f sec., time_each_probe:%.2f sec\n', ...
            time_in_parallel, time_in_serial, time_each_probe);
end
diary off;

