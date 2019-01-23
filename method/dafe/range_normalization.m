function f = range_normalization(f)
% save('./temp/range_normalization.mat', 'F');

% clear
% clc
% load('./temp/range_normalization.mat');

n = size(f,2);
f_max = repmat(max(f),n,1);
f_min = repmat(min(f),n,1);
f = (f-f_min)./(f_max - f_min);







