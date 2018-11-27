function [id_tab, name_tab] = generate_feedback_suggestion(f, V, labeled_gallery_ix, gallery_name_tab, ctrl_para)
% save('./temp/generate_feedback_suggestion.mat', 'f', 'V', 'labeled_gallery_ix', 'gallery_name_tab', 'ctrl_para');

% clear
% clc
% load('./temp/suggest_feedback_list.mat');

K = size(V,2);
fbppr = ctrl_para.exp_para.fbppr;
method = ctrl_para.dafe.feedback_suggestion_method;
gallery_set_num = ctrl_para.dataset.gallery_set_num;
id_tab = zeros(fbppr, K);
name_tab = cell(fbppr,K);

f = (1+f)/2; % normalize f to [0,1] so that it has the same range as V
F = repmat(f,1,K);
rank_F = zeros(size(F));

for k=1:K
    labeled_gallery_ix_k = labeled_gallery_ix;
    labeled_gallery_ix_k(labeled_gallery_ix_k==(gallery_set_num+1)) = [];
    nl = length(labeled_gallery_ix_k);
    if nl>0
        F(labeled_gallery_ix_k,k) = -1*ones(nl,1);
    end
    [~, ix] = sort(F(:,k), 'descend');
    [~, rank_F(:,k)] = sort(ix);
end
    
switch method
    case 'V-only'
        A = V;

    case 'f-only'
        A = F;
        
    case 'V+f'
        A = V + F;
        
    case 'V*f'
        A = V.*F;
        
    case 'V/rank(f)'
        
        A = V./rank_F;
        
    case 'rank(V)/rank(f)'
        
        rank_V = zeros(size(V));
        for k=1:K
            [~, ix] = sort(V(:,k), 'ascend');
            [~, rank_V(:,k)] = sort(ix);
        end

        A = rank_V./rank_F;
end

epsilon = 1e-3;
for k=1:K
    [a, ix] = sort(A(:,k), 'descend');
    
    %% this trick is to kept ix result output by different machines are the same
    if a(1)<=1 && length(find(a>1-epsilon))>fbppr
        rng(K); 
        little_disturbance = rand(316,1)*epsilon;
        [~, ix] = sort(A(:,k)+little_disturbance, 'descend');
    end
%     if length(find(a>1-epsilon))>15
%         ix(1:10)'
%         stop = 1;
%     end
    %% trick end
    
    for i=1:fbppr
        id_tab(i,k) = ix(i);
        name_tab{i,k} = gallery_name_tab{id_tab(i,k)};
    end
end

