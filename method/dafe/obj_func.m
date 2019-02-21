function [J, item_smooth, item_fitting, item_pace, item_sparse] = obj_func(f, v, y, W, model_para)


alpha = diag(model_para.alpha);
beta = model_para.beta;
gamma = model_para.gamma;
n = model_para.node_set_num;

v_tilde = repmat(v,1,n)+repmat(v',n,1);
W_tilde = v_tilde.*W;
P_tilde = diag(sum(W_tilde,2));
Q_tilde = diag(sum(v_tilde,2));
P = diag(sum(W,2));

R = P\P_tilde-sqrt(P)\W_tilde/sqrt(P);

item_smooth = 2*f'*R*f/(n*n);
item_fitting = 2*(f-y)'*(alpha*Q_tilde)*(f-y)/(n*n);
item_pace = ones(1,n)*(beta*v_tilde)*ones(n,1)/(n*n);
item_sparse = gamma*v'*v/n;

J = item_smooth+item_fitting-item_pace+item_sparse;



