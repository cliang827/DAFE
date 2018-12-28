function [f, fval] = solve_f(v, y, W, model_para)
% purpose: 
%   solve f directly without splitting it into ulabeled and labeled parts, 
%   i.e., f = [fu; fl]; (This appears in AAAI'19 draft)
%   
% save('./temp/solve_f.mat', 'v', 'y', 'W', 'model_para', 'ctrl_para');

% clear
% clc
% load('./temp/solve_f.mat');

alpha = model_para.alpha;
n = length(y);

v_tilde = (1-v)*(1-v)';
Q_hat = diag(alpha.*sum(v_tilde,2));

P = diag(sum(W,2));
W_tilde = v_tilde.*W;

% figure
% subplot(1,2,1); imshow(W);
% subplot(1,2,2); imshow(W_tilde);

W_hat = sqrt(P)\W_tilde/sqrt(P);
P_tilde = diag(sum(W_tilde,2));
P_hat = sqrt(P)\P_tilde/sqrt(P);

A = (P_hat - W_hat + Q_hat + 1e-6*eye(n));
b = 2*(Q_hat*y); 


A = (A+A')/2;
[~, p]=chol(A); 
assert(0==p); % assure A is positive-definite. Large V (>=1) values may trigger this error

%% matlab solver
% tic
% epsilon = 1e-6;
H = 2*A;
z = -b;
lb = -1*ones(n,1);
ub = ones(n,1);
options = optimoptions('quadprog','Algorithm','interior-point-convex','Display','off');
[f, fval] = quadprog(H,z,[],[],[],[],lb,ub,[],options);
% time = toc;

% assert(norm(f(nu+1:n) - fl)<epsilon);

