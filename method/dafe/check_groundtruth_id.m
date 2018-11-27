function [id_tab, name_tab] = check_groundtruth_id(f, probe_id, query_times, gallery_name_tab, ...
    history_id_tab, history_name_tab, ctrl_para)
% save('./temp/check_groundtruth_id.mat', 'f', 'probe_id', 'query_times', 'gallery_name_tab', 'history_id_tab', 'history_name_tab', 'ctrl_para');

% clear
% clc
% load('./temp/check_groundtruth_id.mat');

id_tab = history_id_tab{query_times};
name_tab = history_name_tab{query_times};
gt_id = [];

switch ctrl_para.dataset.name
    case 'viper'
        % in single shot re-id dataset, gt_id == probe_id
        [~, temp_ix] = sort(f, 'descend');
        [~, rank_f] = sort(temp_ix);
        gt_rank = rank_f(probe_id);
        
        if gt_rank<=ctrl_para.exp_para.rank_threshold
            gt_id = probe_id;
        end 
end

if ~isempty(gt_id)
    K = ctrl_para.dataset.K;
    fbppr = ctrl_para.exp_para.fbppr;
    gt_num = length(gt_id);
    
    if gt_num>0 && ctrl_para.exp_para.include_groundtruth_flag
       
        % generation: gt_id_tab
        for k=1:K
            % check what id we have already labeled
            id_set = [];
            for qt=1:query_times
                id_set = cat(1, id_set, history_id_tab{qt}(:,k));
            end
            
            % identify gt_id that not labeled
            gt_id_tab = [];
            for i=1:gt_num
                if ~ismember(gt_id(i), id_set)
                    gt_id_tab = cat(1, gt_id_tab, gt_id(i));
                end
            end
            
            % forming final id_tab
            temp_id_tab = cat(1, gt_id_tab, history_id_tab{query_times}(:,k));
            id_tab(:,k) = temp_id_tab(1:fbppr);
            
            % forming final name_tab
            name_tab(:,k) = gallery_name_tab(id_tab(:,k));
        end
    end
end

