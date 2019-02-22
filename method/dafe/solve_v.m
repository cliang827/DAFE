function [v, beta] = solve_v(f, v0, y, W, model_para)
% save('./temp/solve_v.mat', 'f', 'v', 'y', 'W', 'model_para');

% clear 
% clc
% load('./temp/solve_v.mat');


alpha = model_para.alpha;
beta_percentage = model_para.beta_percentage;
gamma = model_para.gamma;
fb_num = model_para.fb_num;
v_sum_constraint = model_para.v_sum_constraint;
epsilon = 1e-20;
labeled_gallery_set = model_para.labeled_gallery_set;
unlabeled_gallery_set = model_para.unlabeled_gallery_set;
n = model_para.node_set_num;

P = diag(sum(W,2));
f_normalized = sqrt(P)\f;
ff = repmat(f_normalized,[1 n])-repmat(f_normalized',[n 1]);
alpha_fY = repmat(alpha.*(f-y).*(f-y), [1 n]) + repmat(alpha'.*(f-y)'.*(f-y)',[n 1]);
L = W.*ff.*ff + alpha_fY;

L_temp = L;
L_temp(labeled_gallery_set,:) = []; 
L_temp(:,labeled_gallery_set) = []; 
% % L_temp = sort(L_temp(:), 'descend');
% % L_temp = L_temp(1:2:end);
% % beta = L_temp(max(1,floor(beta_percentage*length(L_temp))));
L_mean = sort(mean(L_temp));
beta = L_mean(max(1, floor(beta_percentage*n)));
            
L_hat = L-beta;
figure(1); a = L_hat;a(a<0) = 0;a(a>0) = 1; subplot(1,2,1); imshow(a);

b = 2*(sum(L_hat,2));
H = 2*(epsilon+gamma*n)*eye(n);

if v_sum_constraint
    Aeq = zeros(2, n);
    Aeq(1,labeled_gallery_set) = 1;
    Aeq(2,unlabeled_gallery_set) = 1;
    beq = [length(labeled_gallery_ix);length(unlabeled_gallery_set)-fb_num];
else
    Aeq = zeros(1, n);
    Aeq(1,labeled_gallery_set) = 1;
    beq = [length(labeled_gallery_set)];
end

lb = zeros(n,1);
ub = ones(n,1);
options = optimoptions('quadprog','Algorithm','interior-point-convex','Display','off');
v = quadprog(H,b,[],[],Aeq,beq,lb,ub,v0,options);


subplot(1,2,2); plot(v);


