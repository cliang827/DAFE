function [cmc_result ] = result_cmc(reid_score, groundtruth)


[~, ordered] = sort(reid_score, 'ascend');  %ÿһ�дӴ�С����
match = (ordered == groundtruth);

% cmc result
cmc_result = cumsum(sum(match, 2)./size(match, 2));







