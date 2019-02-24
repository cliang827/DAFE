clear
clc



% dataset - GRID
% result_mat_file = './result/GRID/mmap-pc-2019-02-24-02-41-32.mat'; 

% dataset - PRID450S
% result_mat_file = './result/PRID450s/mmap-pc-2019-02-23-15-00-49.mat'; 

% dataset - VIPeR
% result_mat_file = './result/VIPeR/mmap-pc-2019-02-23-10-57-19.mat'; 

% dataset - CUHK01
result_mat_file = './result/CUHK01/mmap-pc-2019-02-23-10-54-50.mat'; 

% dataset - CUHK03labeled
% result_mat_file = './result/CUHK03labeled/mmap-pc-2019-02-23-19-47-21.mat';

% dataset - CUHK03detected
% result_mat_file = './result/CUHK03detected/mmap-pc-2019-02-23-14-55-00.mat';

load(result_mat_file);

% each dataset only run one trial
% trial_set = cell2mat(para_test_set(:,7));
% feedback_id = feedback_id(trial_set==trial_set(1),:);
% para_test_set = para_test_set(trial_set==trial_set(1),:);
% assert(sum(cell2mat(para_test_set(:,6))-[1 2 3]')==0);
clearvars -except feedback_id para_test_set eval_para


[status,cmdout] = system('hostname');
if strcmp(eval_para.machine_type,'mmap-pc') && ...
        strcmp(cmdout(1:4), 'x270')
    eval_para.data_file_dir = strrep(eval_para.data_file_dir,'\','/');
end
load(eval_para.data_file_dir);

para_test_num = size(para_test_set,1);
dataset_list = cell(para_test_num,1);
feedback_id_list = cell(para_test_num,1);
para_list = cell(para_test_num,1);

for i=1:para_test_num
    trail_ix = para_test_set{i,7};
    para_list{i}.probe_feat = probe_feat_set{trail_ix};
    para_list{i}.gallery_feat = gallery_feat_set{trail_ix};
    para_list{i}.method_set = {'mr', 'emr'};
    para_list{i}.tot_query_times = max(cell2mat(para_test_set(:,6)));
    
    % parameters for manifold ranking (mr)
    para_list{i}.mr.k = 10;
    para_list{i}.mr.alpha = 0.99;
    
    % parameters for efficient manifold ranking (emr)
    para_list{i}.emr.alpha = 0.99;
    para_list{i}.emr.p = 10;
    
    
    gallery_set_num = size(para_list{i}.gallery_feat,1);
    dataset_list{i}.node_set_num = gallery_set_num + 1;
    robot_dist = robot_dist_set{trail_ix};
    dataset_list{i}.robot_feedback_score = dist2sim(robot_dist);
    dataset_list{i}.groundtruth_rank = groundtruth_rank_set{trail_ix};
    
    feedback_id_list{i} = feedback_id{i};
end

%% open parallel/serial
if ~isempty(gcp('nocreate'))>0
    delete(gcp('nocreate'))
end

poolobj = parpool(feature('NumCores')); %parpool;
batch_size = poolobj.NumWorkers;   


cmc_score_list = cell(para_test_num, 1);
parfor i=1:para_test_num
    [~, cmc_score_list{i}] = baseline(dataset_list{i}, feedback_id_list{i}, para_list{i});
end

if ~isempty(gcp('nocreate'))>0
    delete(gcp('nocreate'))
end

trial_set = unique(cell2mat(para_test_set(:,7)));
page1_table_list = zeros(2,2,length(trial_set));
page2_table_list = zeros(12,2,length(trial_set));
for i=1:para_test_num 
    cmc_score = cmc_score_list{i};
    trial_ix = find(trial_set==para_test_set{i,7});
    fb_num = para_test_set{i,6};
    page1_table_list(:,:,trial_ix) = cat(1, cmc_score.mr([1 20],1)', cmc_score.emr([1 20],1)');
    page2_table_list((2*fb_num-1):(2*fb_num),:,trial_ix) = cat(1, cmc_score.mr([1 20],2)', cmc_score.emr([1 20],2)');
    page2_table_list((6+2*fb_num-1):(6+2*fb_num),:,trial_ix) = cat(1, cmc_score.mr([1 20],3)', cmc_score.emr([1 20],3)');
end

page1_table = mean(page1_table_list,3);
page1_table = 100*page1_table;

page2_table = mean(page2_table_list,3);
page2_table = 100*page2_table;