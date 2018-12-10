% function [reid_score, auc_score, feedback_id, time_result] = baseline(dataset, feedback_id, para)
% save('./temp/baseline.mat', 'dataset', 'feedback_id', 'para');

clear
clc
load('./temp/baseline.mat');


probe_feat = dataset.probe_feat;
gallery_feat = dataset.gallery_feat;
probe_set_num = dataset.probe_set_num;
gallery_set_num = dataset.gallery_set_num;
node_set_num = dataset.node_set_num;
robot_feedback_score = dataset.robot_feedback_score;

method_set = para.method_set;
tot_query_times = para.tot_query_times;

if ismember('emr', method_set)
    reid_score.emr = zeros(gallery_set_num, probe_set_num, tot_query_times);
end
if ismember('mr', method_set)
    reid_score.mr = zeros(gallery_set_num, probe_set_num, tot_query_times);
end



for i=1:probe_set_num
    feat_data = cat(2, gallery_feat, probe_feat(:,i))';
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
        
        
        %% Sigir'11_emr
        if ismember('emr', method_set)
            f = emr(feat_data, y, para.emr); 
            reid_score.emr(:,i,query_times) = f(1:gallery_set_num);
        end
        
        %% NIPS'03_mr
        if ismember('emr', method_set)
            f = mr(feat_data, y, para.mr); 
            reid_score.emr(:,i,query_times) = f(1:gallery_set_num);
        end
    end
end