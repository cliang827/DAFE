function [reid_score, difficulty_score, auc_score, suggest_feedback_id_tab, time_result] = ...
    dafe(dataset, ctrl_para)
% save('./temp/dafe.mat', 'dataset', 'ctrl_para');

% clear 
% clc
% load('./temp/dafe.mat');


show_progress_flag = ctrl_para.exp.show_progress_flag;
tot_query_times = ctrl_para.exp.tot_query_times;

probe_set_num = dataset.probe_set_num;
gallery_set_num = dataset.gallery_set_num;
node_set_num = dataset.node_set_num;
probe_feat = dataset.probe_feat;
gallery_feat = dataset.gallery_feat;
feat_dim = dataset.feat_dim;
gallery_name_tab = dataset.gallery_name_tab;
robot_feedback_score = dist2sim(dataset.robot_dist);

reid_score_y = zeros(gallery_set_num, probe_set_num, tot_query_times);
reid_score_f = zeros(gallery_set_num, probe_set_num, tot_query_times);
reid_score_f_h1 = zeros(gallery_set_num, probe_set_num, tot_query_times);
reid_score_f_h2 = zeros(gallery_set_num, probe_set_num, tot_query_times);
reid_score_f_h3 = zeros(gallery_set_num, probe_set_num, tot_query_times);
reid_score_f_y_method1 = zeros(gallery_set_num, probe_set_num, tot_query_times);
reid_score_f_y_method2 = zeros(gallery_set_num, probe_set_num, tot_query_times);
difficulty_score = zeros(gallery_set_num, probe_set_num, tot_query_times);
suggest_feedback_id_tab = cell(probe_set_num, tot_query_times);
suggest_feedback_name_tab = cell(probe_set_num, tot_query_times);
query_time_tab = zeros(probe_set_num, tot_query_times);
iter_time_tab = zeros(probe_set_num, tot_query_times);

model_para.alpha = ctrl_para.model.alpha;
beta_percentage = ctrl_para.model.beta_percentage;
model_para.gamma = ctrl_para.model.gamma;
model_para.fb_num = ctrl_para.model.fb_num;
model_para.v_sum_constraint = ctrl_para.exp.v_sum_constraint;
model_para.test_history_flag = ctrl_para.exp.test_history_flag;
model_para.test_y_method_flag = ctrl_para.exp.test_y_method_flag;
model_para.node_set_num = dataset.node_set_num;

v0 = ones(node_set_num,1);


show_progress_step = 10;
for i=1:probe_set_num
    if show_progress_flag && (mod(i,show_progress_step)==0 || i==1)
        nchar = fprintf(1, 'progress=%d/%d (%3.0f%%) ...', i, probe_set_num, 100*i/probe_set_num);
    end
    
    node_feat = cat(1, gallery_feat, probe_feat(i,:));
    feedback_scores = [];
    labeled_gallery_set = [];
    
    for query_times = 1:tot_query_times
        start_time = tic;

        % update feedback_gallery_ix and feedback_scores
        if 1==query_times
            new_feedback_scores = 1;
            new_feedback_gallery_ix = gallery_set_num + 1;
        else
            new_feedback_gallery_ix = suggest_feedback_id_tab{i, query_times-1};
            new_feedback_scores = robot_feedback_score(new_feedback_gallery_ix,i);
        end
        feedback_scores  = cat(1, feedback_scores, new_feedback_scores);
        labeled_gallery_set = cat(1, labeled_gallery_set, new_feedback_gallery_ix);
        unlabeled_gallery_set = setdiff(1:node_set_num, labeled_gallery_set);

        
        %% prepare parameters for DAFE
        % parameter: W
        if 1==query_times
            M = eye(feat_dim);
        else
            pairwise_ix_set = generate_pairwise_set(labeled_gallery_set, feedback_scores);            
            M = similarity_learning(M, pairwise_ix_set, node_feat, 'oasis');
        end
        W = node_feat*M*node_feat';
        W(eye(node_set_num)==1) = 0;
        
        % parameter: f0
        f0 = ones(node_set_num,1); 
        f0(1:gallery_set_num) = range_normalization(W(1:gallery_set_num,end));
        
        % parameter: y0
        y0 = zeros(node_set_num, 1);
        y0(node_set_num) = 1;
        
        % parameter: beta
%         if query_times == 1
            P = diag(sum(W,2));
            f_normalized = sqrt(P)\f0; % eq.(31) in TR17
            ff = repmat(f_normalized,[1 node_set_num])-repmat(f_normalized',[node_set_num 1]);
            smooth_loss = W.*ff.*ff;
            fitting_loss = repmat(ctrl_para.model.alpha.*(f0-y0).*(f0-y0), [1 node_set_num]) + ...
                repmat(ctrl_para.model.alpha'.*(f0-y0)'.*(f0-y0)',[node_set_num 1]);
            total_loss = smooth_loss + fitting_loss;
            total_loss(labeled_gallery_set,:) = []; 
            total_loss(:,labeled_gallery_set) = []; 
            sorted_total_loss = sort(total_loss(:), 'descend');
            sorted_total_loss = sorted_total_loss(1:2:end);
            model_para.beta = sorted_total_loss(max(1,floor(beta_percentage*length(sorted_total_loss))));
%         end

        % others parameters
        model_para.labeled_gallery_set = labeled_gallery_set;
        model_para.unlabeled_gallery_set = unlabeled_gallery_set;
        model_para.y_labeled = y0;
        model_para.y_labeled(labeled_gallery_set) = feedback_scores;
        model_para.y_labeled_v2 = f0;
        model_para.y_labeled_v2(labeled_gallery_set) = feedback_scores;

        [f, v, f_y_method, f_history, iter_times] = solve_fv(f0, v0, y0, W, model_para);
        y = f0(1:gallery_set_num);
%         figure(1); plot(v);

        %% result collection
        reid_score_f_y_method1(:,i,query_times) = f_y_method(1:gallery_set_num,1);
        reid_score_f_y_method2(:,i,query_times) = f_y_method(1:gallery_set_num,2);
        reid_score_f_h1(:,i,query_times) = squeeze(f_history(1:gallery_set_num,1));
        reid_score_f_h2(:,i,query_times) = squeeze(f_history(1:gallery_set_num,2));
        reid_score_f_h3(:,i,query_times) = squeeze(f_history(1:gallery_set_num,3));
        reid_score_f(:,i,query_times) = squeeze(f(1:gallery_set_num));
        reid_score_y(:,i,query_times) = squeeze(y(1:gallery_set_num));
        difficulty_score(:,i,query_times) = v(1:gallery_set_num);

        
        %% generate feedback suggestion
        curr_reid_score = reid_score_f(:, i, query_times);
        curr_difficulty_score = difficulty_score(:, i, query_times);

        [suggest_feedback_id_tab{i,query_times}, suggest_feedback_name_tab{i,query_times}] = ...
            generate_feedback_suggestion(curr_reid_score, curr_difficulty_score, ...
            labeled_gallery_set, gallery_name_tab, ctrl_para);

        if ctrl_para.exp.include_groundtruth_flag
            [suggest_feedback_id_tab{i,query_times}, suggest_feedback_name_tab{i,query_times}] = ...
                check_groundtruth_id(curr_reid_score, i, query_times, gallery_name_tab, ...
                suggest_feedback_id_tab(i,:), suggest_feedback_name_tab(i,:), ctrl_para);
        end
        
        iter_time_tab(i,query_times) = iter_times;
        query_time_tab(i, query_times) = toc(start_time);
    end
 
    if show_progress_flag && mod(i,show_progress_step)==show_progress_step-1
        fprintf(1, repmat('\b', 1, nchar));
    end
end    
time_result.time_by_round = mean(query_time_tab,1);
time_result.time_in_total = sum(query_time_tab(:));

%%
reid_score.y = reid_score_y;
reid_score.f = reid_score_f;

auc_score_y = zeros(1, tot_query_times);
auc_score_f_y_method1 = zeros(1, tot_query_times);
auc_score_f_y_method2 = zeros(1, tot_query_times);
auc_score_f = zeros(1, tot_query_times);
auc_score_f_h1 = zeros(1, tot_query_times);
auc_score_f_h2 = zeros(1, tot_query_times);
auc_score_f_h3 = zeros(1, tot_query_times);
groundtruth_rank = dataset.groundtruth_rank;

for query_times = 1:tot_query_times

    [~, auc_score_f_y_method1(1, query_times)] = ...
        result_evaluation(reid_score_f_y_method1(:,:,query_times), groundtruth_rank);
    
    [~, auc_score_f_y_method2(1, query_times)] = ...
        result_evaluation(reid_score_f_y_method2(:,:,query_times), groundtruth_rank);
    
    [~, auc_score_f_h1(1, query_times)] = ...
        result_evaluation(reid_score_f_h1(:,:,query_times), groundtruth_rank);
    
    [~, auc_score_f_h2(1, query_times)] = ...
        result_evaluation(reid_score_f_h2(:,:,query_times), groundtruth_rank);
    
    [~, auc_score_f_h3(1, query_times)] = ...
        result_evaluation(reid_score_f_h3(:,:,query_times), groundtruth_rank);
    
    [~, auc_score_f(1, query_times)] = ...
        result_evaluation(reid_score_f(:,:,query_times), groundtruth_rank);

    [~, auc_score_y(1, query_times)] = ...
        result_evaluation(reid_score_y(:,:,query_times), groundtruth_rank);
end
auc_score.y = auc_score_y;
auc_score.f_y_method1 = auc_score_f_y_method1;
auc_score.f_y_method2 = auc_score_f_y_method2;
auc_score.f = auc_score_f;
auc_score.f_h1 = auc_score_f_h1;
auc_score.f_h2 = auc_score_f_h2;
auc_score.f_h3 = auc_score_f_h3;

%%
if show_progress_flag
    fprintf(1, repmat('\b', 1, nchar));
    fprintf(1, 'progress=%d/%d (100%%) ... time=%.2f sec.\n', probe_set_num, probe_set_num, time_result.time_in_total);
    
    
    qt = [1:tot_query_times]';
    y = 100*auc_score_y';
    f_y_method1 = 100*auc_score_f_y_method1';
    f_y_method2 = 100*auc_score_f_y_method2';
    f = 100*auc_score_f';
    f_h1 = 100*auc_score_f_h1';
    f_h2 = 100*auc_score_f_h2';
    f_h3 = 100*auc_score_f_h3';
    T = table(qt,y,f_y_method1,f_y_method2,f,f_h1,f_h2,f_h3)
    

    for qt = 1:tot_query_times
        if qt==1, fprintf('\n\t\ty\tf_ym1\tf_ym2\tf\tf_h1\tf_h2\n'); end
        fprintf(1, 'qt=%d,\t auc = [%.2f%%\t%.2f%%\t%.2f%%\t%.2f%%|\t%.2f%%\t%.2f%% ]\n', ...
            qt, 100*auc_score_y(qt), 100*auc_score_f_y_method1(qt),  100*auc_score_f_y_method2(qt), ...
            100*auc_score_f(qt), 100*auc_score_f_h1(qt), 100*auc_score_f_h2(qt));
    end
end
