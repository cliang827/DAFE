function [reid_score, auc_score, feedback_id, time_result] = baseline(dataset, feedback_id, para)
% save('./temp/baseline.mat', 'dataset', 'feedback_id', 'para');

% clear
% clc
% load('./temp/baseline.mat');

node_set_num = dataset.node_set_num;
robot_feedback_score = dataset.robot_feedback_score;

method_set = para.method_set;
tot_query_times = para.tot_query_times;
probe_feat = para.probe_feat;
gallery_feat = para.gallery_feat;
probe_set_num = size(probe_feat,1);
gallery_set_num = size(gallery_feat,1);


if ismember('mr', method_set)
    reid_score.mr = zeros(gallery_set_num, probe_set_num, tot_query_times);
    auc_score.mr = zeros(1, tot_query_times);
%     feedback_id.mr = [];
    time_result.mr = zeros(probe_set_num, tot_query_times);
end
if ismember('emr', method_set)
    rng('default');
    reid_score.emr = zeros(gallery_set_num, probe_set_num, tot_query_times);
    auc_score.emr = zeros(1, tot_query_times);
%     feedback_id.emr = [];
    time_result.emr = zeros(probe_set_num, tot_query_times);
end

for i=1:probe_set_num
    feat_data = cat(1, gallery_feat, probe_feat(i,:));
    y = zeros(node_set_num,1);
        
    for query_times = 1:tot_query_times
        if 1==query_times
            new_feedback_scores = 1;
            new_feedback_gallery_ix = gallery_set_num + 1;
        else
            new_feedback_gallery_ix = feedback_id{i, query_times-1};
            new_feedback_scores = robot_feedback_score(new_feedback_gallery_ix,i);
        end
        y(new_feedback_gallery_ix) = new_feedback_scores;
        
        %% NIPS'03_mr     
        if ismember('mr', method_set)
            t_start = tic;
            f = mr(feat_data, y, para.mr); 
            reid_score.mr(:,i,query_times) = f(1:gallery_set_num);
            time_result.mr(i,query_times) = toc(t_start);
        end
        
        %% Sigir'11_emr
        if ismember('emr', method_set)
            t_start = tic;
            f = emr(feat_data, y, para.emr); 
            reid_score.emr(:,i,query_times) = f(1:gallery_set_num);
            time_result.emr(i,query_times) = toc(t_start);
        end
    end
end

%%result analysis
groundtruth_rank = repmat(1:probe_set_num, gallery_set_num, 1);
if ismember('mr', method_set)
    for query_times = 1:tot_query_times
        [~, auc_score.mr(1, query_times)] = ...
            result_evaluation(reid_score.mr(:,:,query_times), groundtruth_rank);
    end
    time_result.time_by_round.mr = mean(time_result.mr,1);
    time_result.time_in_total.mr = sum(time_result.mr(:));
    time_result = rmfield(time_result, 'mr');
end
if ismember('emr', method_set)
    for query_times = 1:tot_query_times
        [~, auc_score.emr(1, query_times)] = ...
            result_evaluation(reid_score.emr(:,:,query_times), groundtruth_rank);
    end
    time_result.time_by_round.emr = mean(time_result.emr,1);
    time_result.time_in_total.emr = sum(time_result.emr(:));
    time_result.emr = [];
    time_result = rmfield(time_result, 'emr');
end
stop = 1;
