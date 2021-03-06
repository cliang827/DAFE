function [id_tab, name_tab] = generate_feedback_suggestion(f, v, labeled_gallery_ix, gallery_name_tab, ctrl_para)
% save('./temp/generate_feedback_suggestion.mat', 'f', 'v', 'labeled_gallery_ix', 'gallery_name_tab', 'ctrl_para');

% clear
% clc
% load('./temp/suggest_feedback_list.mat');

fb_num = ctrl_para.model.fb_num;
fb_method = ctrl_para.model.fb_method;
gallery_set_num = length(f);
id_tab = zeros(fb_num,1);
name_tab = cell(fb_num,1);

f = (1+f)/2; % normalize f to [0,1] so that it has the same range as v
labeled_gallery_ix(labeled_gallery_ix==(gallery_set_num+1)) = [];
nl = length(labeled_gallery_ix);
if nl>0
    f(labeled_gallery_ix) = -1*ones(nl,1);
end
[~, ix] = sort(f, 'descend');
[~, rank_f] = sort(ix);

feedback_score = zeros(size(f));
switch fb_method
    case 'v-only'
        feedback_score = v;

    case 'f-only'
        feedback_score = f;
        
    case 'v+f'
        feedback_score = v + f;
        
    case 'v*f'
        feedback_score = v.*f;
        
    case 'v/rank(f)'
        feedback_score = v./rank_f;
        
    case 'rank(v)/rank(f)'
        [~, ix] = sort(v, 'ascend');
        [~, rank_v] = sort(ix);
        feedback_score = rank_v./rank_f;
end

[a, ix] = sort(feedback_score, 'descend');
epsilon = 1e-3;
%% this trick is to kept ix result output by different machines are the same
if a(1)<=1 && length(find(a>1-epsilon))>fb_num
    rng(K); 
    little_disturbance = rand(316,1)*epsilon;
    [~, ix] = sort(feedback_score+little_disturbance, 'descend');
end
%% trick end

for i=1:fb_num
    id_tab(i) = ix(i);
    name_tab{i} = gallery_name_tab{id_tab(i)};
end


