%% prepare dataset and directories
curr_dataset.name = 'viper';
curr_dataset.source = 'cvpr16_gog_xqda';

ctrl_para.dir_info.result_dir = ['.' slash 'result' slash];
if ~exist(ctrl_para.dir_info.result_dir, 'dir'), mkdir(ctrl_para.dir_info.result_dir); end

ctrl_para.dir_info.temp_dir = ['.' slash 'temp' slash];
if ~exist(ctrl_para.dir_info.temp_dir, 'dir'), mkdir(ctrl_para.dir_info.temp_dir); end

ctrl_para.dir_info.data_dir = ['.' slash 'data' slash curr_dataset.name slash];
if ~exist(ctrl_para.dir_info.data_dir, 'dir'), error('no data file!'); end

ctrl_para.dir_info.data_file = [ctrl_para.dir_info.data_dir slash curr_dataset.source '.mat'];
data_file = load(ctrl_para.dir_info.data_file);

version_str = cellstr(datetime('now','Format','y-MM-d-HH-mm-ss'));
ctrl_para.dir_info.result_file = sprintf('%s%s-%s.mat', ...
    ctrl_para.dir_info.result_dir, machine_type, version_str{1});
ctrl_para.dir_info.log_file = [ctrl_para.dir_info.result_file(1:end-4),'.txt'];

ctrl_para.dir_info.method_dir = ['.' slash 'method' slash];
 

%% set search ranges of model and experiment parameters
ctrl_para.exp.fb_method_set = {'rank(v)/rank(f)'}; 
ctrl_para.exp.alpha_set = 10.^(0); %10.^(-1:0.5:1);
ctrl_para.exp.beta_percentage_set = 0.05; %[0.05 0.1 0.5 1]; 
ctrl_para.exp.gamma_set = 0; 
ctrl_para.exp.delta_set = 0; %[0.01 0.5 0.99];
if debug_flag
    ctrl_para.exp.fb_num_set = [1 5];
    ctrl_para.exp.trial_set = [1 2];
    ctrl_para.exp.show_progress_flag = true;
    ctrl_para.exp.show_figure_flag = true;
else
    ctrl_para.exp.fb_num_set = 1:10;
    ctrl_para.exp.trial_set = 1:10;
    ctrl_para.exp.show_progress_flag = false;
    ctrl_para.exp.show_figure_flag = false;
end
ctrl_para.exp.tot_query_times = 3; 
ctrl_para.exp.show_table_flag = true;

ctrl_para.exp.v_sum_constraint = false;
ctrl_para.exp.include_groundtruth_flag = false;
ctrl_para.exp.rank_threshold = 20;
ctrl_para.exp.machine_type = machine_type;
ctrl_para.exp.run_mode = run_mode;
ctrl_para.exp.batch_size = batch_size;


ctrl_para.model.tau = 0.1;                               % used for PCM'14 initialization
ctrl_para.model.p = 0.9;
ctrl_para.model.regu_method = 'cvpr07_spectral_matting'; %'cvpr07_spectral_matting';

%% construct parameter grid for searching the optimal configuration
method_num = length(ctrl_para.exp.fb_method_set);
alpha_num = length(ctrl_para.exp.alpha_set);
beta_num = length(ctrl_para.exp.beta_percentage_set);
gamma_num = length(ctrl_para.exp.gamma_set);
delta_num = length(ctrl_para.exp.delta_set);
fb_num = length(ctrl_para.exp.fb_num_set);
trial_num = length(ctrl_para.exp.trial_set);

para_test_set = cell(method_num*alpha_num*beta_num*gamma_num*delta_num*fb_num*trial_num, 7);
para_test_num = 0;
for i = 1:method_num
    for j=1:alpha_num
        for k=1:beta_num
            for l=1:gamma_num
                for m=1:delta_num
                    for n=1:fb_num
                        for t=1:trial_num
                            para_test_num = para_test_num + 1;

                            para_test_set{para_test_num,1} = ctrl_para.exp.fb_method_set{i};
                            para_test_set{para_test_num,2} = ctrl_para.exp.alpha_set(j);
                            para_test_set{para_test_num,3} = ctrl_para.exp.beta_percentage_set(k);
                            para_test_set{para_test_num,4} = ctrl_para.exp.gamma_set(l);
                            para_test_set{para_test_num,5} = ctrl_para.exp.delta_set(m);
                            para_test_set{para_test_num,6} = ctrl_para.exp.fb_num_set(n);
                            para_test_set{para_test_num,7} = ctrl_para.exp.trial_set(t);
                        end
                    end
                end
            end 
        end
    end
end


% prepare parameters for parallel implementation
ctrl_para_set = cell(1, para_test_num);
dataset_set = cell(1,para_test_num);
for i=1:para_test_num
    ctrl_para.model.fb_method = para_test_set{i,1};
    ctrl_para.model.alpha = para_test_set{i,2};
    ctrl_para.model.beta_percentage = para_test_set{i,3};
    ctrl_para.model.gamma = para_test_set{i,4};
    ctrl_para.model.delta = para_test_set{i,5};
    ctrl_para.model.fb_num = para_test_set{i,6};
    ctrl_para.exp.trial = para_test_set{i,7};
    ctrl_para_set{i} = ctrl_para;
    
    t = ctrl_para.exp.trial;
    probe_ix_set = data_file.testinds_set{t}(data_file.testcamIDs_set{t}==1);
    gallery_ix_set = data_file.testinds_set{t}(data_file.testcamIDs_set{t}==2);
    curr_dataset.probe_feat = data_file.feature_all(probe_ix_set,:)';
    curr_dataset.gallery_feat = data_file.feature_all(gallery_ix_set,:)';
    dist_id_converter = data_file.dist_id_converter;
    new_gallery_ix_set = dist_id_converter(gallery_ix_set);
    new_probe_ix_set = dist_id_converter(probe_ix_set);
    curr_dataset.g2g_dist = data_file.g2g_dist(new_gallery_ix_set, new_gallery_ix_set);
    curr_dataset.g2p_dist = data_file.g2p_dist(new_gallery_ix_set, new_probe_ix_set);
    [curr_dataset.gallery_set_num, curr_dataset.probe_set_num] = size(curr_dataset.g2p_dist);
    if debug_flag, curr_dataset.probe_set_num = 10; end        
    curr_dataset.node_set_num = curr_dataset.gallery_set_num+1;
    curr_dataset.gallery_name_tab = data_file.allimagenames(gallery_ix_set);
    curr_dataset.robot_feedback_score = data_file.feedback_score_set{t};
    dataset_set{i} = curr_dataset;
end


%% prepare parameters for result analysis
eval_para.result_file = ctrl_para.dir_info.result_file;
eval_para.show_figure_flag = ctrl_para.exp.show_figure_flag;
eval_para.show_table_flag = ctrl_para.exp.show_table_flag;
eval_para.trial_num = trial_num;
eval_para.v_sum_constraint = ctrl_para.exp.v_sum_constraint;
eval_para.tot_query_times = ctrl_para.exp.tot_query_times;
eval_para.probe_set_num = curr_dataset.probe_set_num;
eval_para.gallery_set_num = curr_dataset.gallery_set_num;
eval_para.machine_type = machine_type;
eval_para.fb_num_set = ctrl_para.exp.fb_num_set;
eval_para.data_file = ctrl_para.dir_info.data_file;