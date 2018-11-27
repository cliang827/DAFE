function [cmc_result, auc_result, rank_result, hist_result] = result_evaluation(reid_score, groundtruth)
% save('.\temp\result_evaluation.mat', 'reid_score', 'groundtruth');

% clc 
% clear all
% close all
% load('.\temp\result_evaluation.mat');


[~, ordered] = sort(reid_score, 'descend'); 
match = (ordered == groundtruth);

% cmc result
cmc_result = cumsum(sum(match, 2)./size(match, 2));

% auc_result
auc_result = cmc2auc(cmc_result);
% auc_result = 0.5*(2*sum(cmc_result) - cmc_result(1) - cmc_result(end))/(size(reid_score, 1)-1);

% hist_result
[rank_result, ~] = find(match);
hist_result = hist(rank_result,1:length(rank_result));










