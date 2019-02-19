function [cmc_result ] = result_cmc(reid_score, groundtruth, method)


[~, ordered] = sort(reid_score, method);  
match = (ordered == groundtruth);

% cmc result
cmc_result = cumsum(sum(match, 2)./size(match, 2));







