function v = solve_v(f, v0, y, W, model_para)
% save('./temp/solve_v.mat', 'f', 'v', 'y', 'W', 'model_para');

% clear 
% clc
% load('./temp/solve_v.mat');


alpha = model_para.alpha;
beta = model_para.beta;
gamma = model_para.gamma;
epsilon = 1e-9;
labeled_gallery_ix = model_para.labeled_gallery_set;
n = model_para.node_set_num;

P = diag(sum(W,2));
f_normalized = sqrt(P)\f;
ff = repmat(f_normalized,[1 n])-repmat(f_normalized',[n 1]);
alpha_fY = repmat(alpha.*(f-y).*(f-y), [1 n]) + repmat(alpha'.*(f-y)'.*(f-y)',[n 1]);
L = W.*ff.*ff + alpha_fY;

L_hat = L-beta;

b = 2*(sum(L_hat,2));
H = 2*(epsilon+gamma)*n*eye(n);

Aeq = zeros(1, n);
Aeq(1,labeled_gallery_ix) = 1;
beq = length(labeled_gallery_ix);

lb = zeros(n,1);
ub = ones(n,1);
options = optimoptions('quadprog','Algorithm','interior-point-convex','Display','off');
v = quadprog(H,b,[],[],Aeq,beq,lb,ub,v0,options);





