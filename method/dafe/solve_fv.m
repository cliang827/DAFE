function [f, v, f_y_method, f_history, iter_times] = solve_fv(f, v, y, W, model_para)
% save('./temp/solve_fv.mat', 'f', 'v', 'y', 'W', 'model_para');


% clear
% clc
% load('./temp/solve_fv.mat');


iter_times = 0;
epsilon_J = 1e-6;
outer_loop_max_iter_times = 3;

node_set_num = model_para.node_set_num;
% v_one = ones(node_set_num,1);
y_labeled = model_para.y_labeled;
y_labeled_v2 = model_para.y_labeled_v2;

J_val = zeros(2, outer_loop_max_iter_times);
f_y_method = zeros(node_set_num, 2);
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
    f = solve_f(f, v, y, W, model_para);
    J_val(2, iter_times) = obj_func(f, v, y, W, model_para);
    f_history(:,iter_times) = f;

    % test different y methods
    if iter_times==1 && model_para.test_y_method_flag
        f_y_method(:,1) = solve_f([], v, y_labeled, W, model_para);

        f_y_method(:,2) = solve_f([], v, y_labeled_v2, W, model_para);
    end

    if model_para.test_history_flag
        % test multi-round re-ranking results
        if iter_times>outer_loop_max_iter_times
            break;
        end
    else
        if iter_times>=outer_loop_max_iter_times || ...
                abs(J_val(1, iter_times)-J_val(2, iter_times))<epsilon_J
            break;
        end
    end
end
    




