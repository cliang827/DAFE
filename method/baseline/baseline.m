function [reid_score, auc_score, feedback_id, time_result] = baseline(dataset, feedback_id, para)
% save('./temp/baseline.mat', 'dataset', 'feedback_id', 'para');

% clear
% clc
% load('./temp/baseline.mat');

node_set_num = dataset.node_set_num;
robot_feedback_score = dist2sim(dataset.robot_dist);


method_set = para.method_set;
method_set_icmr17 = para.method_set_icmr17;
method_set_other = para.method_set_other;
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

%% icmr17 7 Algorithms
% 'CPRR','RLRECOM','RKGRAPH','CORGRAPH','RLSIM','RECKNNGRAPH','CONTEXTRR'
for method_idx =1:length(method_set_icmr17)
    cur_icmr_method = method_set_icmr17{method_idx};
    
    for i=1:probe_set_num
        
        feat_tmp = cat(1, probe_feat(i,:), gallery_feat);
        dist_tmp = pdist2(feat_tmp, feat_tmp, 'euclidean'); 
        
        for query_times = 1:tot_query_times
            
            dataset_size = size(feat_tmp, 1);
            [score, icmr_time] = icmr17(dist_tmp, cur_icmr_method, dataset_size);
            reid_score.icmr17(:,i,query_times, method_idx) = score;
            % update distance for later reranking
            dist_tmp(1, 2:end) = -log(score);
            dist_tmp(2:end, 1) = -log(score);  
            time_result.icmr17(i, query_times, method_idx) = icmr_time;
            
        end
    end 
end

%% krnn
for method_idx =1:length(method_set_other)
    feat_data = cat(1, probe_feat, gallery_feat);
    dist_tmp = pdist2(feat_data, feat_data, 'euclidean');
    for query_times = 1:tot_query_times

        t_start = tic;
        dist_reranking = re_ranking(dist_tmp, probe_set_num, para.krnn.k1, ...
            para.krnn.k2, para.krnn.lambda);
        time_result.krnn(query_times) = toc(t_start);

        dist_tmp = dist_reranking; 
        reid_score.krnn(:,:,query_times) = exp(-dist_reranking(1:probe_set_num,probe_set_num+1:end)');  
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

%% icmr17 eval
for method_idx =1:length(method_set_icmr17)
    for query_times = 1:tot_query_times
        [~, auc_score.icmr17(method_idx, query_times)] = ...
            result_evaluation(reid_score.icmr17(:,:,query_times, method_idx), groundtruth_rank);
    end
    time_result.time_by_round.icmr17(method_idx, :) = mean(time_result.icmr17(:,:,method_idx),1);
    time_result.time_in_total.icmr17(method_idx, :) = sum(sum(time_result.icmr17(:,:, method_idx)),2); 
end
if isfield(time_result, 'icmr17'), time_result = rmfield(time_result,'icmr17');end


%% krnn eval
for method_idx =1:length(method_set_other)
    for query_times = 1:tot_query_times
        [~, auc_score.krnn(1, query_times)] = ...
            result_evaluation(reid_score.krnn(:,:,query_times), groundtruth_rank);
    end
    time_result.time_by_round.krnn = time_result.krnn./probe_set_num;
    time_result.time_in_total.krnn = sum(time_result.krnn(:));
end
if isfield(time_result, 'krnn'), time_result = rmfield(time_result, 'krnn');end

stop = 1;
