function [Eq_12, Eq_15, Eq_18, Eq_21, Eq_22] = ...
    Eq_func2(v, L_hat, gamma, L_hat_plus, L_hat_minus, d_hat, vt, p, regu_method, epsilon)
save('./temp/Eq_func2.mat', 'v', 'L_hat', 'gamma', 'L_hat_plus', 'L_hat_minus', 'd_hat', 'vt', 'p', 'regu_method', 'epsilon');

% clear
% clc
% load('./temp/Eq_func2.mat');

% Eq.(12)
n = length(v);


% %% test that whether g(v) is a convex function
% g_v  =  v'*L_hat_minus*v/(n*n) -gamma*regu(v,p,0,regu_method)/n;
% g_vt = vt'*L_hat_minus*vt/(n*n)-gamma*regu(vt,p,0,regu_method)/n;
% for lambda = 0:0.001:1
%     
%     if lambda==0.069
%         stop = 1;
%     end
%     
%     v_lambda_vt = lambda*v + (1-lambda)*vt;
%     g_v_lambda_vt = v_lambda_vt'*L_hat_minus*v_lambda_vt/(n*n)-gamma*regu(v_lambda_vt,p,0,regu_method)/n;
%     
%     assert(g_v_lambda_vt-(lambda*g_v + (1-lambda)*g_vt)<=epsilon);
% end
% %% end of test



Eq_12_term_1 = (1-v)'*L_hat*(1-v)/(n*n);
Eq_12_term_2 = gamma*regu(v,p,0,regu_method)/n;
Eq_12 = Eq_12_term_1 + Eq_12_term_2;

% Eq.(15)
fv = v'*L_hat_plus*v/(n*n)-2*v'*d_hat/(n*n);
gv = v'*L_hat_minus*v/(n*n)-gamma*regu(v,p,0,regu_method)/n;
Eq_15 = fv - gv;

L_hat = L_hat_plus-L_hat_minus;
const_eq_13 = sum(L_hat(:))/(n*n);
assert(abs(Eq_12-(Eq_15 + const_eq_13))<epsilon);



% Eq.(20)
partial_norm_vt = regu(vt,p,1,regu_method);


% Eq.(19)
partial_g_vt = 2*L_hat_minus*vt/(n*n) - gamma*partial_norm_vt/n;

% Eq.(18)
g_vt = vt'*L_hat_minus*vt/(n*n)-gamma*regu(vt,p,0,regu_method)/n;
Eq_18 = g_vt + partial_g_vt'*(v-vt);
% note:
%   if error happens in the following assert, it shows that g(v) is not convex
%   to fix this problem, make a big enough initilization of L_hat_plus so
%   that L_hat_minux = L_hat_plus - L_hat is PSD enough to make g(v) is
%   convex
% assert(Eq_18-gv<=epsilon);

% Eq.(21)
const_eq_21 = vt'*L_hat_minus*vt/(n*n)-gamma*regu(vt,p,0,regu_method)/n - ...
    vt'*(2*L_hat_minus*vt/(n*n) - gamma*partial_norm_vt/n);
Eq_21_term_1 = v'*(2*L_hat_minus*vt/(n*n) - gamma*partial_norm_vt/n);
Eq_21_term_2 = const_eq_21;
Eq_21 = Eq_21_term_1 + Eq_21_term_2;
assert(abs(Eq_18-Eq_21)<epsilon);


% Eq.(22)
Eq_22_term_1 = fv;
Eq_22_term_2 = Eq_21;
Eq_22 = Eq_22_term_1 - Eq_22_term_2;












