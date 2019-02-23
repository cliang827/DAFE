function [paras_filtered, cmc_f_filtered, cmc_y_filtered, auc_y_filtered, auc_f_y_method1_filtered, auc_f_y_method2_filtered, ...
    auc_f_filtered, auc_f_h1_filtered, auc_f_h2_filtered, auc_f_h3_filtered, ...
    auc_mr_filtered, auc_emr_filtered, ...
    v_max_filtered, v_mean_filtered, v_std_filtered] = ...
    para_filter(paras, cmc_f, cmc_y, auc_y, auc_f_y_method1, auc_f_y_method2, ...
    auc_f, auc_f_h1, auc_f_h2, auc_f_h3, auc_mr, auc_emr, ...
    v_max, v_mean, v_std, filter)
% save('./temp/para_filter.mat', 'paras', 'cmc_f', 'cmc_y', 'auc_y', 'auc_f_y_method1', 'auc_f_y_method2', ...
%     'auc_f', 'auc_f_h1', 'auc_f_h2', 'auc_f_h3', 'auc_mr', 'auc_emr', ...
%     'v_max', 'v_mean', 'v_std', 'filter');

% clear
% clc
% load('./temp/para_filter.mat');

% if strcmp(filter.column.name_set, 'feedback_method')
%     assert(size(paras,1)==1);
%     paras_filtered = paras;
%     auc_f_y_method_filtered = auc_f_y_method;
%     auc_f_filtered = auc_f;
%     auc_y_filtered = auc_y;
%     v_max_filtered = v_max;
%     v_mean_filtered = v_mean;
%     v_std_filtered = v_std;    
%     return
% end

%% row filter
invalid_para_ix = zeros(size(paras,1),1);
if isfield(filter, 'row') 
    if isfield(filter.row, 'feedback_method_set')
        invalid_feedback_method_set = setdiff(unique(paras(:,1)),filter.row.feedback_method_set);
        for i=1:length(invalid_feedback_method_set)
            invalid_para_ix(paras(:,1)==invalid_feedback_method_set(i))=1;    
        end
    end
    if isfield(filter.row, 'alpha_set')
        invalid_alpha_set = setdiff(unique(paras(:,2)),filter.row.alpha_set);
        for i=1:length(invalid_alpha_set)
            invalid_para_ix(paras(:,2)==invalid_alpha_set(i))=1;    
        end
    end
    if isfield(filter.row, 'beta_set')
        invalid_beta_set = setdiff(unique(paras(:,3)),filter.row.beta_set);
        for i=1:length(invalid_beta_set)
            invalid_para_ix(paras(:,3)==invalid_beta_set(i))=1;
        end
    end
    if isfield(filter.row, 'gamma_set')
        invalid_gamma_set = setdiff(unique(paras(:,4)),filter.row.gamma_set);
        for i=1:length(invalid_gamma_set)
            invalid_para_ix(paras(:,4)==invalid_gamma_set(i))=1;
        end
    end
    if isfield(filter.row, 'delta_set')
        invalid_delta_set = setdiff(unique(paras(:,5)),filter.row.delta_set);
        for i=1:length(invalid_delta_set)
            invalid_para_ix(paras(:,5)==invalid_delta_set(i))=1;
        end
    end
    if isfield(filter.row, 'fb_num_set')
        invalid_fbppr_set = setdiff(unique(paras(:,6)),filter.row.fbppr_set);
        for i=1:length(invalid_fbppr_set)
            invalid_para_ix(paras(:,6)==invalid_fbppr_set(i))=1;
        end
    end
end

paras(invalid_para_ix==1,:) = [];
cmc_f(invalid_para_ix==1,:,:) = [];
cmc_y(invalid_para_ix==1,:,:) = [];
auc_y(invalid_para_ix==1,:) = [];
auc_f_y_method1(invalid_para_ix==1,:) = [];
auc_f_y_method2(invalid_para_ix==1,:) = [];
auc_f(invalid_para_ix==1,:) = [];
auc_f_h1(invalid_para_ix==1,:) = [];
auc_f_h2(invalid_para_ix==1,:) = [];
auc_f_h3(invalid_para_ix==1,:) = [];
auc_mr(invalid_para_ix==1,:) = [];
auc_emr(invalid_para_ix==1,:) = [];
v_max(invalid_para_ix==1,:) = [];
v_mean(invalid_para_ix==1,:) = [];
v_std(invalid_para_ix==1,:) = [];

%% column filter
filter_column_ix = [];
if ~isempty(strfind(filter.column.name_set, 'feedback_method'))
    filter_column_ix = cat(2, filter_column_ix, 1);
end
if ~isempty(strfind(filter.column.name_set, 'alpha'))
    filter_column_ix = cat(2, filter_column_ix, 2);
end
if ~isempty(strfind(filter.column.name_set, 'beta'))
    filter_column_ix = cat(2, filter_column_ix, 3);
end
if ~isempty(strfind(filter.column.name_set, 'gamma'))
    filter_column_ix = cat(2, filter_column_ix, 4);
end
if ~isempty(strfind(filter.column.name_set, 'delta'))
    filter_column_ix = cat(2, filter_column_ix, 5);
end
if ~isempty(strfind(filter.column.name_set, 'fb_num'))
    filter_column_ix = cat(2, filter_column_ix, 6);
end

   

if ~isempty(filter_column_ix)
    paras_set = unique(paras(:,filter_column_ix),'stable', 'rows');
    n = size(paras_set,1);
    paras_filtered = zeros(n, size(paras,2));
    
    cmc_f_filtered = zeros(n, size(cmc_f,2),size(cmc_f,3));
    cmc_y_filtered = zeros(n, size(cmc_y,2),size(cmc_y,3));
    auc_y_filtered = zeros(n, size(auc_y,2));
    auc_f_y_method1_filtered = zeros(n, size(auc_f_y_method1,2));
    auc_f_y_method2_filtered = zeros(n, size(auc_f_y_method2,2));
    auc_f_filtered = zeros(n, size(auc_f,2));
    auc_f_h1_filtered = zeros(n, size(auc_f_h1,2));
    auc_f_h2_filtered = zeros(n, size(auc_f_h2,2));
    auc_f_h3_filtered = zeros(n, size(auc_f_h3,2));
    auc_mr_filtered = zeros(n, size(auc_mr,2));
    auc_emr_filtered = zeros(n, size(auc_emr,2));
    v_max_filtered = zeros(n, size(v_max,2));
    v_mean_filtered = zeros(n, size(v_mean,2));
    v_std_filtered = zeros(n, size(v_std,2));

    for i=1:n
        temp = paras(:,filter_column_ix)-repmat(paras_set(i,:), size(paras,1), 1);
        temp2 = mean(abs(temp),2);
        filter_ix_set = find(temp2==0);

        paras_filtered(i,:) = mean(paras(filter_ix_set,:),1);
        cmc_f_filtered(i,:,:) = mean(cmc_f(filter_ix_set, :,:),1);
        cmc_y_filtered(i,:,:) = mean(cmc_y(filter_ix_set, :,:),1);
        auc_y_filtered(i,:) = mean(auc_y(filter_ix_set, :),1);
        auc_f_y_method1_filtered(i,:) = mean(auc_f_y_method1(filter_ix_set, :),1);
        auc_f_y_method2_filtered(i,:) = mean(auc_f_y_method2(filter_ix_set, :),1);
        auc_f_filtered(i,:) = mean(auc_f(filter_ix_set, :),1);
        auc_f_h1_filtered(i,:) = mean(auc_f_h1(filter_ix_set, :),1);
        auc_f_h2_filtered(i,:) = mean(auc_f_h2(filter_ix_set, :),1);
        auc_f_h3_filtered(i,:) = mean(auc_f_h3(filter_ix_set, :),1);
        auc_mr_filtered(i,:) = mean(auc_mr(filter_ix_set, :),1);
        auc_emr_filtered(i,:) = mean(auc_emr(filter_ix_set, :),1);
        v_max_filtered(i,:) = mean(v_max(filter_ix_set, :),1);
        v_mean_filtered(i,:) = mean(v_mean(filter_ix_set,:),1);
        v_std_filtered(i,:) = mean(v_std(filter_ix_set,:),1);    
    end
else
    n = size(paras,1);
    if n==1
        return;
    end
    
    error('needs revision!');
    
%     paras_filtered = mean(paras,1);
%     auc_y_filtered = mean(auc_y,1);
%     auc_f_y_method1_filtered = mean(auc_f_y_method1,1);
%     auc_f_y_method2_filtered = mean(auc_f_y_method2,1);
%     auc_f_filtered = mean(auc_f,1);
%     auc_f_h1_filtered = mean(auc_f_h1,1);
%     auc_f_h2_filtered = mean(auc_f_h2,1);
%     auc_f_h3_filtered = mean(auc_f_h3,1);
%     auc_mr_filtered = mean(auc_mr,1);
%     auc_emr_filtered = mean(auc_emr,1);
%     v_max_filtered = mean(v_max,1);
%     v_mean_filtered = mean(v_mean,1);
%     v_std_filtered = mean(v_std,1);    
end





