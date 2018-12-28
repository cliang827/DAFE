function [score] = mr(fea,y0,opts)
% [score] = emr(data,y0,opts): Manifold Ranking
% Input:
%       - data: the data matrix of size nSmp x nFea, where each row is a sample
%               point
%       - y0:   the initial query vector, e.g., query item =1 and the other all 0; 
%       
%       opts: options for this algorithm
%           - p: the number of landmarks
nSmp = size(fea,1);
W = constructW(fea, opts);
D_mhalf = full(sum(W,2).^-.5);
D_mhalf = spdiags(D_mhalf,0,nSmp,nSmp);
S = D_mhalf*W*D_mhalf;
alpha = opts.alpha;
S = speye(nSmp)-alpha*S; 
score = S\y0;