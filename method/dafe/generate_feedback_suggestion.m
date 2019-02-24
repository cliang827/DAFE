function [id_tab, name_tab] = generate_feedback_suggestion(f, v, labeled_gallery_ix, gallery_name_tab, ctrl_para)
% save('./temp/generate_feedback_suggestion.mat', 'f', 'v', 'labeled_gallery_ix', 'gallery_name_tab', 'ctrl_para');

% clear
% clc
% load('./temp/generate_feedback_suggestion.mat');


rank_threshold = ctrl_para.exp.rank_threshold;
fb_num = ctrl_para.model.fb_num;
fb_method = ctrl_para.model.fb_method;
gallery_set_num = length(f);
id_tab = zeros(fb_num,1);
name_tab = cell(fb_num,1);


labeled_gallery_ix(labeled_gallery_ix==(gallery_set_num+1)) = [];
nl = length(labeled_gallery_ix);

if nl>0
    f(labeled_gallery_ix) = 0;
    v(labeled_gallery_ix) = 1;
end
[~, ix] = sort(f, 'descend');
[~, rank_f] = sort(ix);

feedback_score = zeros(size(f));
epsilon = 1e-6;

switch fb_method
    case 'rank(v)-in-top-k'
        ix_set = ix(1:rank_threshold);
        v_value = v(ix_set);
        [~, ix_temp] = sort(v_value, 'descend');
        [~, rank_v] = sort(ix_temp);
        feedback_score(ix_set) = rank_v;

    case 'v'
        feedback_score = v;
        
    case '1-v'
        feedback_score = 1-v;

    case 'f'
        feedback_score = f;
        
    case 'f+v'
        feedback_score = f + v;
        
    case 'f-v'
        feedback_score = f - v;
        
    case {'v*f', 'f*v'}
        feedback_score = v.*f;
        
    case 'v/rank(f)'
        feedback_score = v./rank_f;
        
    case 'rank(v)/rank(f)'
        [~, ix] = sort(v, 'ascend');
        [~, rank_v] = sort(ix);
        feedback_score = rank_v./(rank_f-10*epsilon);
        
    case 'rank(1-v)/rank(f)'
        [~, ix] = sort(v, 'descend');
        [~, rank_v] = sort(ix);
        feedback_score = rank_v./(rank_f-10*epsilon);
        
    case 'rank(1-v)'
        [~, ix] = sort(v, 'descend');
        [~, rank_v] = sort(ix);
        feedback_score = rank_v;
        
    case 'rank(f)'
        feedback_score = 1./rank_f;
        
    case 'rand'
        feedback_score = rand(length(v),1);
end

[a, ix] = sort(feedback_score, 'descend');
%% this trick is to kept ix result output by different machines are the same
if a(1)-a(1+fb_num)<epsilon
    disp('no discrimination of current feedback score!');
    rng('default'); 
    trunc_num = sum(a>(a(1)-epsilon));
    sorted_ix = setdiff(sort(ix(1:trunc_num)),labeled_gallery_ix);
    feedback_score(sorted_ix) = feedback_score(sorted_ix) + ...
        rand(length(sorted_ix),1)*epsilon;
    [~, ix] = sort(feedback_score, 'descend');
end
%% trick end

for i=1:fb_num
    id_tab(i) = ix(i);
    name_tab{i} = gallery_name_tab{id_tab(i)};
end


