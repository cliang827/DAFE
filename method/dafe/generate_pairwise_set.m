function [pairwise_set] = generate_pairwise_set(feedback_gallery_ix, feedback_scores)
% save('./temp/generate_pairwise_set.mat', 'feedback_gallery_ix', 'feedback_scores');

% clear
% clc
% load('./temp/update_pairwise_set.mat');

n = length(feedback_gallery_ix);
compare_set = nchoosek(1:n,2);
compare_set_num = size(compare_set,1);
pairwise_set = zeros(compare_set_num,3);
for set=1:compare_set_num
    ix1 = feedback_gallery_ix(compare_set(set,1));
    ix2 = feedback_gallery_ix(compare_set(set,2));
    if feedback_scores(compare_set(set,1))>=feedback_scores(compare_set(set,2))
        score_diff = feedback_scores(compare_set(set,1))-feedback_scores(compare_set(set,2));
        pairwise_set(set,:) = [ix1, ix2, score_diff];
    else
        score_diff = feedback_scores(compare_set(set,2))-feedback_scores(compare_set(set,1));
        pairwise_set(set,:) = [ix2, ix1, score_diff];
    end    
end



