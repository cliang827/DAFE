% function test_dafe_parallel(machine_type)

clear
close all
clc

ctrl_para.exp_para.machine_type = 'x270';


ctrl_para.dataset.probe_set_num = 316;   % 60
ctrl_para.exp_para.tot_query_times = 3; 
ctrl_para.exp_para.v_sum_constraint = 0;

ctrl_para.exp_para.fb_sugg_method_set = {'rank(V)/rank(f)'}; 
ctrl_para.exp_para.alpha_set = 10.^(0); %10.^(-1:0.5:1);
ctrl_para.exp_para.beta_percentage_set = 0.05; %[0.05 0.1 0.5 1]; 
ctrl_para.exp_para.gamma_set = 0; 
ctrl_para.exp_para.delta_set = 0; %[0.01 0.5 0.99];
ctrl_para.exp_para.fbppr_set = 3; %[1,2,3,5];


ctrl_para.dataset.name = 'viper';
ctrl_para.dataset.gallery_set_num = 316;
ctrl_para.dataset.node_set_num = ctrl_para.dataset.gallery_set_num+1;
ctrl_para.dataset.K = 1; 
ctrl_para.exp_para.user_type = 'robot';
ctrl_para.exp_para.feature_type = 'cvpr16_gog'; % cvpr16_gog, pcm14_riro
ctrl_para.exp_para.robot_name = 'cvpr16_gog_xqda_train';
ctrl_para.exp_para.data_file_name = 'cvpr16_gog_xqda';
ctrl_para.exp_para.trial_num = 10;   % 5

ctrl_para.exp_para.include_groundtruth_flag = false;
ctrl_para.exp_para.rank_threshold = 20;
ctrl_para.exp_para.rep_date_str = '2018-10-26';
ctrl_para.exp_para.inc_date_str = '27-Oct-2018';
ctrl_para.exp_para.show_progress_flag = false;
ctrl_para.exp_para.show_figure_flag = false;
ctrl_para.exp_para.show_table_flag = true;

ctrl_para.dafe.tau = 0.1;           % used for PCM'14 initialization
ctrl_para.dafe.p = 0.9;
ctrl_para.dafe.regu_method = 'cvpr07_spectral_matting'; %'cvpr07_spectral_matting';

switch ctrl_para.exp_para.machine_type
    case {'whuhpc'}
        setenv('TZ','Asia/Shanghai');
        cd('/data/liangchao/ScaDAFE');
        ctrl_para.dir_info.slash ='/';
%         ctrl_para.dir_info.vlfeat_dir = '/data/liangchao/vlfeat-0.9.21/';
%         run([ctrl_para.dir_info.vlfeat_dir 'toolbox/vl_setup']);
        
    case {'laptop', 'cliang-X270', 'x270'}
        cd('/home/cliang/work/code/test/reid/DAFE');
        ctrl_para.dir_info.slash ='/';
%         ctrl_para.dir_info.vlfeat_dir = '/home/cliang/work/code/toolbox/vlfeat-0.9.21/';
%         run([ctrl_para.dir_info.vlfeat_dir 'toolbox/vl_setup']);

    case {'pc', 'mmap', 'mmap-pc', 'cliang-mmap-pc'}
        cd('D:\work\code\test\ScaDAFE');
        ctrl_para.dir_info.slash ='\';
%         ctrl_para.dir_info.vlfeat_dir = 'D:\work\code\test\vlfeat-0.9.21\';
%         run([ctrl_para.dir_info.vlfeat_dir 'toolbox\vl_setup']);
end

slash = ctrl_para.dir_info.slash;
ctrl_para.dir_info.result_dir = ['.' slash 'result' slash];
ctrl_para.dir_info.method_dir = ['.' slash 'method' slash];
ctrl_para.dir_info.temp_dir = ['.' slash 'temp' slash];
ctrl_para.dir_info.data_dir = ['.' slash 'data' slash ctrl_para.dataset.name slash];
ctrl_para.dir_info.data_file = [ctrl_para.dir_info.data_dir slash ctrl_para.exp_para.data_file_name '.mat'];
data_file = load(ctrl_para.dir_info.data_file);
addpath(genpath(fullfile(ctrl_para.dir_info.method_dir))); 


version_str = cellstr(datetime('now','Format','y-MM-d-HH-mm-ss'));
result_file = sprintf('%s%s-%s.mat', ctrl_para.dir_info.result_dir, ...
    ctrl_para.exp_para.machine_type, version_str{1});
diary([result_file(1:end-4),'.txt']); diary on;

fprintf(1, 'dataset=%s, probe_set_num=%d, trial_num=%d, tot_query_times=%d\n', ...
    ctrl_para.dataset.name, ctrl_para.dataset.probe_set_num, ...
    ctrl_para.exp_para.trial_num, ctrl_para.exp_para.tot_query_times);

fprintf(1, 'include_groundtruth_flag=%d, v_sum_constraint_flag=%d\n', ...
    ctrl_para.exp_para.include_groundtruth_flag, ctrl_para.exp_para.v_sum_constraint);

fprintf(1, 'feature_type=''%s'', robot_name=''%s''\n', ...
    ctrl_para.exp_para.feature_type, ctrl_para.exp_para.robot_name);

method_num = length(ctrl_para.exp_para.fb_sugg_method_set);
alpha_num = length(ctrl_para.exp_para.alpha_set);
beta_num = length(ctrl_para.exp_para.beta_percentage_set);
gamma_num = length(ctrl_para.exp_para.gamma_set);
delta_num = length(ctrl_para.exp_para.delta_set);
fbppr_num = length(ctrl_para.exp_para.fbppr_set);
trial_num = ctrl_para.exp_para.trial_num;

para_test_set = cell(method_num*alpha_num*beta_num*gamma_num*delta_num*fbppr_num*trial_num, 7);
para_test_num = 0;
for i = 1:method_num
    for j=1:alpha_num
        for k=1:beta_num
            for l=1:gamma_num
                for m=1:delta_num
                    for n=1:fbppr_num
                        for t=1:trial_num
                            para_test_num = para_test_num + 1;

                            para_test_set{para_test_num,1} = ctrl_para.exp_para.fb_sugg_method_set{i};
                            para_test_set{para_test_num,2} = ctrl_para.exp_para.alpha_set(j);
                            para_test_set{para_test_num,3} = ctrl_para.exp_para.beta_percentage_set(k);
                            para_test_set{para_test_num,4} = ctrl_para.exp_para.gamma_set(l);
                            para_test_set{para_test_num,5} = ctrl_para.exp_para.delta_set(m);
                            para_test_set{para_test_num,6} = ctrl_para.exp_para.fbppr_set(n);
                            para_test_set{para_test_num,7} = t;
                        end
                    end
                end
            end 
        end
    end
end

run_mode_flag = 'parallelization';
if strcmp(run_mode_flag, 'parallelization')==1 
    if ~isempty(gcp('nocreate'))>0
        delete(gcp('nocreate'))
    end

    poolobj = parpool(feature('NumCores')); %parpool;
    batch_size = poolobj.NumWorkers;   
else
    batch_size = 1;
end

ctrl_para_set = cell(1, para_test_num);
dataset_set = cell(1,para_test_num);
for i=1:para_test_num
    ctrl_para.dafe.feedback_suggestion_method = para_test_set{i,1};
    ctrl_para.dafe.alpha = para_test_set{i,2};
    ctrl_para.dafe.beta_percentage = para_test_set{i,3};
    ctrl_para.dafe.gamma = para_test_set{i,4};
    ctrl_para.dafe.delta = para_test_set{i,5};
    ctrl_para.exp_para.fbppr = para_test_set{i,6};
    ctrl_para.exp_para.t = para_test_set{i,7};
    ctrl_para.exp_para.run_mode_flag = run_mode_flag;
    ctrl_para.exp_para.batch_size = batch_size;
    ctrl_para_set{i} = ctrl_para;
    
    t = ctrl_para.exp_para.t;
    probe_ix_set = data_file.testinds_set{t}(data_file.testcamIDs_set{t}==1);
    gallery_ix_set = data_file.testinds_set{t}(data_file.testcamIDs_set{t}==2);
    dist_id_converter = data_file.dist_id_converter;
    new_gallery_ix_set = dist_id_converter(gallery_ix_set);
    new_probe_ix_set = dist_id_converter(probe_ix_set);
    dataset.g2g_dist = data_file.g2g_dist(new_gallery_ix_set, new_gallery_ix_set);
    dataset.g2p_dist = data_file.g2p_dist(new_gallery_ix_set, new_probe_ix_set);
    dataset.gallery_name_tab = data_file.allimagenames(gallery_ix_set);
    dataset.robot_feedback_score = data_file.feedback_score_set{t};
    dataset_set{i} = dataset;
end


reid_score = cell(para_test_num,1);
difficulty_score = cell(para_test_num,1);
feedback_id = cell(para_test_num,1);
time_result = cell(para_test_num,1);
t_start =tic;
parfor i=1:para_test_num
        [reid_score{i}, difficulty_score{i}, feedback_id{i}, time_result{i}] = ...
            dafe(dataset_set{i}, ctrl_para_set{i});
end
if ~isempty(gcp('nocreate'))>0
    delete(gcp('nocreate'))
end
time_in_total_parallel = toc(t_start);

save(result_file, 'para_test_set', 'ctrl_para', 'reid_score', 'difficulty_score', 'feedback_id', 'time_result');
[time_in_total_serial, time_each_probe] = analyze_parameters(result_file);
fprintf(1, '\n time = [%.0f, %0.f, %.2f] sec.\n', time_in_total_parallel, time_in_total_serial, time_each_probe);
diary off;

