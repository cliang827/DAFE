function [J, item_smooth, item_fitting, item_pace, item_sparse] = obj_func(f, v, y, W, model_para)


alpha = diag(model_para.alpha);
beta = model_para.beta;
gamma = model_para.gamma;
p = model_para.p;
regu_method = model_para.regu_method;

v_tilde = (1-v)*(1-v)';
W_tilde = v_tilde.*W;
P_tilde = diag(sum(W_tilde,2));
Q_tilde = diag(sum(v_tilde,2));
P = diag(sum(W,2));

R = P\P_tilde-sqrt(P)\W_tilde/sqrt(P);
m = length(f);

item_smooth = 2*f'*R*f/(m*m);
item_fitting = 2*(f-y)'*(alpha*Q_tilde)*(f-y)/(m*m);
item_pace = ones(1,m)*(beta*v_tilde)*ones(m,1)/(m*m);
item_sparse = gamma*regu(v,p,0,regu_method)/m;

J = item_smooth+item_fitting-item_pace+item_sparse;



