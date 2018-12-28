function [v, fval_rec, L] = solve_v(f, v, y, W, model_para)
% save('./temp/solve_v.mat', 'f', 'v', 'y', 'W', 'model_para');

% clear
% clc
% load('./temp/solve_for_v3.mat');

v_sum_constraint_flag = model_para.v_sum_constraint_flag;
alpha = model_para.alpha;
beta = model_para.beta;
gamma = model_para.gamma;
delta = model_para.delta;
p = model_para.p;
regu_method = model_para.regu_method;
expected_feedback_num = model_para.expected_feedback_num;
labeled_gallery_ix = model_para.labeled_gallery_set;
unlabeled_gallery_ix = model_para.unlabeled_gallery_set;
n = model_para.node_set_num;


epsilon = 1e-6; 
inner_loop_max_iter_times = 5;

fval_rec = zeros(11+expected_feedback_num, inner_loop_max_iter_times);
% 1-5th dimensions:     Eq.(12), Eq.(15), Eq.(18), Eq.(21), Eq.(22)
% 6th dimension:        Obj(vt) = vt'*A*vt-vt'*b
% 7th dimension:        Obj(v) = v'*A*v-v'*b
% 8th dimension:        norm(v-vt)
% 9th dimension:        sparseness
% 10th dimension:       Obj value without regu 
% 11th dimension:       regu value
% 12th - end dimension:  idx of expected_feedback_num largest v
    
P = diag(sum(W,2));
f_normalized = sqrt(P)\f;
ff = repmat(f_normalized,[1 n])-repmat(f_normalized',[n 1]);
alpha_fY = repmat(alpha.*(f-y).*(f-y), [1 n]) + repmat(alpha'.*(f-y)'.*(f-y)',[n 1]);
L = W.*ff.*ff + alpha_fY;


% L_value = sort(L(:), 'descend');
% beta = L_value(floor(length(L_value)*0.05));

L_hat = L-beta;

% temp = L_hat;
% temp(temp>0) = 0;
% temp(temp<0) = 255;
% imshow(temp);




% L_hat(eye(n)==1) = 0;

% assert(0==norm(L_hat-L_hat'));
d_hat = sum(L_hat,2);
eig_value = abs(eig(L_hat));
lambda = max(eig_value)+epsilon;
% lambda = n*n;

L_hat_plus = lambda*eye(n);
L_hat_minus = L_hat_plus-L_hat;
[~, pp]=chol(L_hat_minus); 
assert(0==pp); 

A = L_hat_plus; 
vt = v;
sum_of_v = expected_feedback_num; %sum(v);
iter_times = 0;
while 1
    iter_times = iter_times+1;
    
%     if iter_times == 14
%         stop = 1;
%     end
    
    if norm(vt)<epsilon
        partial_vt = zeros(n,1);
    else
        partial_vt = regu(vt,p,1,regu_method);
    end
    
    b = 2*d_hat+2*L_hat_minus*vt-n*gamma*partial_vt;
    if v_sum_constraint_flag
        c = zeros(2, n);
        c(1,labeled_gallery_ix) = 1;
        c(2,unlabeled_gallery_ix) = 1;
        e = [0;sum_of_v];
    else
        c(1,labeled_gallery_ix) = 1;
        e = 0;
    end
    
    
    
    if iter_times==1
        v = vt;
    else
       %% solve for v
%         tic
        H = 2*A;
        z = -b;
        Aeq = c;
        beq = e;
        lb(unlabeled_gallery_ix) = delta;
        lb(labeled_gallery_ix) = 0;
%         ub = ones(n,1)-epsilon;
%         lb = zeros(n,1);
        ub = ones(n,1);
        
        options = optimoptions('quadprog','Algorithm','interior-point-convex','Display','off');
        v = quadprog(H,z,[],[],Aeq,beq,lb,ub,[],options);
%         time = time + toc;

        if max(v)<1e-10
            stop = 1;
        end
    end
    sparseness = (sqrt(n)-norm(v,1)/norm(v,2))/(sqrt(n)-1);
    assert(~isempty(v));
    
    [~, ix] = sort(v, 'descend');
    [Eq_12, Eq_15, Eq_18, Eq_21, Eq_22, ] = ...
        Eq_func2(v, L_hat, gamma, L_hat_plus, L_hat_minus, d_hat, vt, p, regu_method, epsilon);
    fval_rec(:,iter_times) = [  Eq_12; ... 
                                Eq_15; ... 
                                Eq_18; ...
                                Eq_21; ...
                                Eq_22; ...
                                vt'*A*vt-vt'*b; ...
                                v'*A*v-v'*b; ...
                                norm(v-vt)/n; ...
                                sparseness; ...
                                (1-v)'*L_hat*(1-v)/(n*n); ...
                                gamma*regu(v,p,0,regu_method)/n; ...
                                ix(1:expected_feedback_num)];

    if (iter_times>1 && fval_rec(8,iter_times)<epsilon) || ...
            iter_times>=inner_loop_max_iter_times
        break;
    else
        vt = v;
    end
end
fval_rec = fval_rec(:,1:iter_times);

% figure
% plot(1:iter_times, fval_rec(1, 1:iter_times), 'r-'); hold on;
% plot(1:iter_times, fval_rec(2, 1:iter_times), 'g--'); hold on;
% plot(1:iter_times, fval_rec(3, 1:iter_times), 'b.-'); hold on;
% name = {'Eq(12)', 'Eq(15)', 'Eq(18)'};
% legend(name,'Location','northeast');
% [vec_sparsity(v(1:316)), max(v(1:316)), min(v(1:316)), median(v(1:316))]
% figure
% histogram(v(1:316));



