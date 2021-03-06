function [f, v, f_mr, f_history] = solve_fv(f, v, y, W, model_para)
% purpose: solve f and v
%      (1) do f-v iteration since qt=1;
% save('./temp/solve_fv.mat', 'f', 'v', 'y', 'W', 'model_para');


% clear
% clc
% load('./temp/solve_fv.mat');


iter_times = 0;
epsilon_J = 1e-6;
outer_loop_max_iter_times = 3;

labeled_gallery_set = model_para.labeled_gallery_set;
node_set_num = model_para.node_set_num;

J_val = zeros(2, outer_loop_max_iter_times);
f_mr = zeros(node_set_num, 2);
f_history = zeros(node_set_num,outer_loop_max_iter_times+1);
v_history = zeros(node_set_num,outer_loop_max_iter_times+1);

f_history(:,end) = f;
v_history(:,end) = v;

while 1
    iter_times = iter_times + 1;

    %v-step
    v = solve_v(f, v, y, W, model_para);
    J_val(1, iter_times) = obj_func(f, v, y, W, model_para);
    v_history(:,iter_times) = v;

    %f-step
    f = solve_f(v, y, W, model_para);
    J_val(2, iter_times) = obj_func(f, v, y, W, model_para);
    f_history(:,iter_times) = f;
    
    % manifold ranking
    if iter_times==1
        v_zero = zeros(node_set_num,1);
        f_mr(:,1) = solve_f(v_zero, y, W, model_para);

        y_labeled = zeros(size(y));
        y_labeled(labeled_gallery_set) = y(labeled_gallery_set);
        f_mr(:,2) = solve_f(v_zero, y_labeled, W, model_para);
    end

    break;
    if iter_times>=outer_loop_max_iter_times || ...
            abs(J_val(1, iter_times)-J_val(2, iter_times))<epsilon_J            
        break;
    end
end
    




