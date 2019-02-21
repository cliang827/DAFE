function [v, fval_rec, L] = solve_v(f, v, y, W, model_para)
% save('./temp/solve_v.mat', 'f', 'v', 'y', 'W', 'model_para');

% clear 
% clc
% load('./temp/solve_v.mat');
debug_flag = 0;

v_sum_constraint_flag = model_para.v_sum_constraint_flag;
alpha = model_para.alpha;
% beta = model_para.beta;
beta = 0;
p = model_para.p;
regu_method = model_para.regu_method;
expected_feedback_num = model_para.expected_feedback_num;
labeled_gallery_ix = model_para.labeled_gallery_set;
unlabeled_gallery_ix = model_para.unlabeled_gallery_set;
n = model_para.node_set_num;

epsilon = 1e-6; 
gamma_cav = model_para.gamma;
gamma_cvx = 0;

if debug_flag
    gamma_cav = 0;
    gamma_cvx = 0;
end

P = diag(sum(W,2));
f_normalized = sqrt(P)\f;
ff = repmat(f_normalized,[1 n])-repmat(f_normalized',[n 1]);
alpha_fY = repmat(alpha.*(f-y).*(f-y), [1 n]) + repmat(alpha'.*(f-y)'.*(f-y)',[n 1]);
L = W.*ff.*ff + alpha_fY;

% if debug_flag
%     L_temp = L;
%     L_temp(labeled_gallery_ix,:) = [];
%     L_temp(:,labeled_gallery_ix) = [];
%     sorted_L_temp = sort(L_temp(:), 'descend');
%     sorted_L_temp = sorted_L_temp(1:2:end);
%     beta = sorted_L_temp(max(1,floor(beta_percentage*length(sorted_L_temp))));
% end


L_hat = L-beta;

if debug_flag
    a = L_hat;
    a(a<0) = 0;
    a(a>0) = 1;
    figure(1); subplot(2,1,1); imshow(a);
end

L_hat_hat = L_hat + gamma_cvx*n*eye(n);

if debug_flag
    [~, pp]=chol(L_hat_hat); 
    if pp==0
        H = 2*L_hat_hat;
        b = -2*(sum(L_hat,2));
        if v_sum_constraint_flag
            Aeq = zeros(2, n);
            Aeq(1,labeled_gallery_ix) = 1;
            Aeq(2,unlabeled_gallery_ix) = 1;
            beq = [0;expected_feedback_num];
        else
            Aeq = zeros(1, n);
            Aeq(1,labeled_gallery_ix) = 1;
            beq = 0;
        end
        lb = zeros(n,1);
        ub = ones(n,1);
        options = optimoptions('quadprog','Algorithm','interior-point-convex','Display','off');
        H=(H+H')/2;
        v_cvx = quadprog(H,b,[],[],Aeq,beq,lb,ub,[],options);
        subplot(2,1,2); plot(v);
    end
end

d_hat = sum(L_hat,2);
eig_value = abs(eig(L_hat_hat));
lambda = max(eig_value);

L_hat_plus = (lambda+epsilon)*eye(n);
L_hat_minus = L_hat_plus-L_hat_hat;
[~, pp]=chol(L_hat_minus); 
assert(0==pp); 


inner_loop_max_iter_times = 5;
fval_rec = zeros(3, inner_loop_max_iter_times);
v_history = zeros(n, inner_loop_max_iter_times);

iter_times = 0;
vt = v;
while 1
    iter_times = iter_times+1;    
    if norm(vt)<epsilon
        partial_vt = zeros(n,1);
    else
        partial_vt = regu(vt,p,1,regu_method);
    end
    
    H = 2*L_hat_plus;
    b = -1*(2*d_hat+2*L_hat_minus*vt+n*gamma_cav*partial_vt);
    if v_sum_constraint_flag
        Aeq = zeros(2, n);
        Aeq(1,labeled_gallery_ix) = 1;
        Aeq(2,unlabeled_gallery_ix) = 1;
        beq = [0;expected_feedback_num];
    else
        Aeq = zeros(1, n);
        Aeq(1,labeled_gallery_ix) = 1;
        beq = 0;
    end
    lb = zeros(n,1);
    ub = ones(n,1);
    options = optimoptions('quadprog','Algorithm','interior-point-convex','Display','off');
    v = quadprog(H,b,[],[],Aeq,beq,lb,ub,[],options);
    assert(~isempty(v));
    
    v_history(:,iter_times) = v;
    fval_rec(1, iter_times) = norm(v-vt)/n;
    fval_rec(2, iter_times) = v'*L_hat_plus*v+f'*v;
    fval_rec(3, iter_times) = (1-v)'*L_hat*(1-v)+...
        gamma_cvx*n*norm(v,2)^2-gamma_cav*n*regu(v,p,0,regu_method);
    
    if (fval_rec(1,iter_times)<epsilon) || iter_times>=inner_loop_max_iter_times
        break;
    else
        vt = v;
    end    
end
fval_rec = fval_rec(:,1:iter_times);

if debug_flag
    figure(2); plot(v);
end


