function [time_total, time_each_probe] = analyze_parameters(result_mat_file, eval_para)
save('./temp/analyze_parameters.mat', 'result_mat_file', 'eval_para');

% close all
% clear
% clc
% load('./temp/analyze_parameters.mat');

load(result_mat_file);

show_figure_flag = eval_para.show_figure_flag;
show_table_flag = eval_para.show_table_flag;
trial_num = eval_para.trial_num;
% v_sum_constraint = eval_para.v_sum_constraint;
tot_query_times = eval_para.tot_query_times;
probe_set_num = eval_para.probe_set_num;
gallery_set_num = eval_para.gallery_set_num;
% machine_type = eval_para.machine_type; 
fb_num_set = eval_para.fb_num_set;
data_file = eval_para.data_file;


feedback_method = para_test_set{1,1};
paras = cell2mat(para_test_set(:,2:6));
paras_num = size(paras,1);
if trial_num>1
    paras = unique(paras,'rows');
    paras_num = paras_num/trial_num;
    auc_score = reshape(auc_score, trial_num, paras_num)';
    difficulty_score = reshape(difficulty_score, trial_num, paras_num)';
    feedback_id = reshape(feedback_id, trial_num, paras_num)';
    time_result = reshape(time_result, trial_num, paras_num)';
end
    

auc_f_mr1_temp = zeros(paras_num,tot_query_times,trial_num);
auc_f_mr2_temp = zeros(paras_num,tot_query_times,trial_num);
auc_f_h1_temp = zeros(paras_num,tot_query_times,trial_num);
auc_f_h2_temp = zeros(paras_num,tot_query_times,trial_num);
auc_f_h3_temp = zeros(paras_num,tot_query_times,trial_num);
auc_f_temp = zeros(paras_num,tot_query_times,trial_num);
auc_y_temp = zeros(paras_num,tot_query_times,trial_num);
v_max_temp = zeros(paras_num,tot_query_times,trial_num);
v_mean_temp = zeros(paras_num,tot_query_times,trial_num);
v_std_temp = zeros(paras_num,tot_query_times,trial_num);
time_round_temp = zeros(paras_num,tot_query_times,trial_num);

auc_f_mr1 = zeros(paras_num,tot_query_times);
auc_f_mr2 = zeros(paras_num,tot_query_times);
auc_f_h1 = zeros(paras_num,tot_query_times);
auc_f_h2 = zeros(paras_num,tot_query_times);
auc_f_h3 = zeros(paras_num,tot_query_times);
auc_f = zeros(paras_num,tot_query_times);
auc_y = zeros(paras_num,tot_query_times);
v_max = zeros(paras_num,tot_query_times);
v_mean = zeros(paras_num,tot_query_times);
v_std = zeros(paras_num,tot_query_times);
time_total = 0;
time_round = zeros(paras_num,tot_query_times);
for i=1:paras_num
    for t=1:trial_num
        time_total = time_total+time_result{i,t}.time_in_total;
        time_round_temp(i,:,t) = time_result{i,t}.time_by_round;
        
        auc_y_temp(i,:,t) = auc_score{i,t}.y;
        auc_f_mr1_temp(i,:,t) = auc_score{i,t}.f_mr1;
        auc_f_mr2_temp(i,:,t) = auc_score{i,t}.f_mr2;
        auc_f_temp(i,:,t) = auc_score{i,t}.f;
        auc_f_h1_temp(i,:,t) = auc_score{i,t}.f_h1;
        auc_f_h2_temp(i,:,t) = auc_score{i,t}.f_h2;
        auc_f_h3_temp(i,:,t) = auc_score{i,t}.f_h3;

        v = squeeze(mean(difficulty_score{i,t},2)); %mean for all probes
        v_max_temp(i,:,t) = max(v);
        v_mean_temp(i,:,t) = mean(v);
        v_std_temp(i,:,t) = std(v);
    end
    
    auc_y(i,:) = mean(auc_y_temp(i,:,:),3);
    auc_f_mr1(i,:) = mean(auc_f_mr1_temp(i,:,:),3);
    auc_f_mr2(i,:) = mean(auc_f_mr2_temp(i,:,:),3);
    auc_f(i,:) = mean(auc_f_temp(i,:,:),3);
    auc_f_h1(i,:) = mean(auc_f_h1_temp(i,:,:),3);
    auc_f_h2(i,:) = mean(auc_f_h2_temp(i,:,:),3);
    auc_f_h3(i,:) = mean(auc_f_h3_temp(i,:,:),3);
    v_max(i,:) = mean(v_max_temp(i,:,:),3);
    v_mean(i,:) = mean(v_mean_temp(i,:,:),3);
    v_std(i,:) = mean(v_std_temp(i,:,:),3);
    time_round(i,:) = mean(time_round_temp(i,:,:),3);
    
end
time_each_probe = time_total/(paras_num*tot_query_times*probe_set_num*trial_num);

%% filter
filter.column.name_set = 'fb_num'; %'feedback_method'; %'alpha-beta-gamma', 'fb_num';
filter.row.alpha_set = 10.^(0);
filter.row.beta_set = 0.05;
[paras, auc_y, auc_f_mr1, auc_f_mr2, auc_f, auc_f_h1, auc_f_h2, auc_f_h3, v_max, v_mean, v_std] = ...
    para_filter(paras, auc_y, auc_f_mr1, auc_f_mr2, auc_f, auc_f_h1, auc_f_h2, auc_f_h3, v_max, v_mean, v_std, filter);

alpha = log10(paras(:,1));
alpha(alpha==-Inf) = -10;
beta = paras(:,2);
gamma = log10(paras(:,3));
gamma(gamma==-Inf) = -10;
delta = paras(:,4);
fbppr = paras(:,5);


%%
result = cat(2, paras, auc_f(:,1), auc_f(:,tot_query_times));
[~, ix_star] = max(result(:,end));
fprintf(1, '\n\nbest parameter set: log(alpha)=%.1f, beta-percentage=%.2f%%, log(gamma)=%.2f, delta=%.2f, fbppr=%.2f\n', ...
    alpha(ix_star), 100*beta(ix_star), gamma(ix_star), delta(ix_star), fbppr(ix_star));
fprintf(1, 'best result: auc1=%.2f%%, auc%d=%.2f%%, better_num=%d/%d\n', ...
    100*auc_f(ix_star,1), tot_query_times, 100*auc_f(ix_star,tot_query_times), ...
    sum(auc_f(:,tot_query_times)>auc_f(:,1)), size(paras,1));

switch filter.column.name_set
    case 'alpha-beta-gamma'
        % figure 1: auc vs. alpha+beta+gamma (3D) and alpha-beta/alpha-gamma/beta-gamma (2D)
        draw_fig_auc_vs_alpha_beta_gamma;

        % figure 2: alpha/beta/gamma vs. auc(f) and mean(v)/std(v)
        draw_fig_alpha_beta_gamma_vs_auc_v;

        % figure 3-x: fix alpha, inspect beta and gamma (qt=1,2,3) 
        draw_fig_v_vs_beta_gamma_qt;
        
    case 'delta'
        hfig = figure;
        delta_set = ctrl_para.exp.delta_set;
        plot(delta_set, mean(auc_y_kmean,2), 'ko--'); hold on;
        plot(delta_set, mean(auc_f_mr,2), 'bs-.'); hold on;
        plot(delta_set, mean(auc_f,2), 'r*-'); hold on; grid on;
        set(gca, 'xtick', delta_set);
        set(gca,'xticklabel',num2cell(delta_set));
        xlabel('delta');
        ylabel('auc');
        legend({'y-kmean', 'f-mr', 'f-kmean'}, 'location', 'southeast');
        saveas(hfig, [result_mat_file(1:end-4), '-delta.fig']);
        if ~show_figure_flag, close(hfig); end

    case 'fb_num'
        hfig = figure;
        ymin = floor(100*min([auc_y(:);auc_f_mr1(:);auc_f_mr2(:);auc_f(:)]))/100;
        ymax = ceil(100*max([auc_y(:);auc_f_mr1(:);auc_f_mr2(:);auc_f(:)]))/100;
        for qt=1:tot_query_times
            subplot(1,tot_query_times, qt);
            plot(fb_num_set, auc_y(:,qt), 'ko--'); hold on;
            plot(fb_num_set, auc_f_mr1(:,qt), 'bs-.'); hold on;
            plot(fb_num_set, auc_f_mr1(:,qt), 'gs-.'); hold on;
            plot(fb_num_set, auc_f(:,qt), 'r*-'); hold on; grid on;
            axis([-Inf Inf ymin ymax]);
            set(gca, 'xtick', fb_num_set);
            set(gca,'xticklabel',num2cell(fb_num_set));
            xlabel('fb\_num');
            ylabel('auc');
            legend({'y', 'f-mr1', 'f-mr2', 'f'}, 'location', 'northwest');
            title(sprintf('qt=%d',qt));
        end
        saveas(hfig, [result_mat_file(1:end-4), '-fb_num.fig']);
        if ~show_figure_flag, close(hfig); end
        
    case 'feedback_method'
        %% method comparison - table version
        auc_result = cat(1, auc_y, auc_f_mr1, auc_f_mr2, auc_f, auc_f_h1, auc_f_h2, auc_f_h3);
        qt = [1:tot_query_times]';
        y = 100*auc_result(1,:)';
        f_mr1 = 100*auc_result(2,:)';
        f_mr2 = 100*auc_result(3,:)';
        f = 100*auc_result(4,:)';
        f_h1 = 100*auc_result(5,:)';
        f_h2 = 100*auc_result(6,:)';
        f_h3 = 100*auc_result(7,:)';
        method_cmp_tab = table(qt,y,f_mr1,f_mr2,f,f_h1,f_h2,f_h3);
        if show_table_flag, method_cmp_tab, end
        
        %% method comparison - figure version
        hfig = figure;
        plot(auc_y, 'ko--'); hold on;
        plot(auc_f_mr1, 'bs-.'); hold on;
        plot(auc_f_mr2, 'g^-.'); hold on;
        plot(auc_f, 'r*-'); hold on; grid on;
        set(gca, 'xtick', 1:tot_query_times);
        set(gca,'xticklabel',{'qt=1', 'qt=2', 'qt=3'});
        xlabel('query times');
        ylabel('auc');
        legend({'y', 'f-mr1', 'f-mr2', 'f'}, 'location', 'southeast');
        saveas(hfig, [result_mat_file(1:end-4), '-method_cmp.fig']);
        if ~show_figure_flag, close(hfig); end
        
        %% feedback details - table version
        if size(paras,1)==1 && trial_num==1
            % LC: for condition of only one set of parameters!
            rank_detail = cell(probe_set_num, 2*tot_query_times+1);
            tab_column_name = cell(1, 2*tot_query_times+1);
            groundtruth_rank = repmat(1:probe_set_num, gallery_set_num, 1);
            for qt=1:tot_query_times
                if qt==1
                    [~, ~, temp] = result_evaluation(reid_score{1,1}.y(:,:,qt), groundtruth_rank); 
                    rank_detail(:,1) = num2cell(temp);

                    [~, ~, temp] = result_evaluation(reid_score{1,1}.f(:,:,qt), groundtruth_rank);
                    rank_detail(:,2) = num2cell(temp);

                    tab_column_name(1,1:2) = {'init', 'qt1'};
                else
                    rank_detail(:,2*qt-1) = feedback_id{1,1}(:,qt-1);

                    [~, ~, temp] = result_evaluation(reid_score{1,1}.f(:,:,qt), groundtruth_rank);
                    rank_detail(:,2*qt) = num2cell(temp);
                    tab_column_name(1,2*qt-1:2*qt) = {sprintf('feedback_ix_rank_score_qt%d',qt-1), sprintf('qt%d',qt)};
                end
            end

            rank_detail(:,end) = num2cell(cell2mat(rank_detail(:,1))-cell2mat(rank_detail(:,2*tot_query_times)));
            tab_column_name{1,end} = 'improvement';

            tab_row_name = cell(probe_set_num,1);
            data_file = load(data_file);
            robot_feedback_score = data_file.feedback_score_set{1};
            for i=1:probe_set_num
                tab_row_name{i} = num2str(i);
                for qt=2:tot_query_times
                    feedback_ix = rank_detail{i,2*qt-1}';
                    [~, ~, feedback_rank] = result_evaluation(repmat(reid_score{1,1}.f(:,i,qt-1),1,fbppr), ...
                        repmat(feedback_ix, gallery_set_num,1));
                    feedback_score = robot_feedback_score(feedback_ix,i);
                    temp = [];
                    for j=1:fbppr
                        if feedback_score(j)>0
                            temp = sprintf('%s: %3d/%d/+%.2f ',temp, feedback_ix(j), feedback_rank(j), feedback_score(j));
                        elseif feedback_score(j)<=0
                            temp = sprintf('%s: %3d/%d/-%.2f ',temp, feedback_ix(j), feedback_rank(j), -1*feedback_score(j));
                        end
                    end
                    rank_detail{i,2*qt-1} = sprintf('[%s]', temp(2:end));
                end
            end
            rank_detail_tab = cell2table(rank_detail,'VariableNames',tab_column_name, 'RowNames', tab_row_name);
            if show_table_flag, rank_detail_tab, end

            %% feedback details - figure version
            hfig = figure;
            suptitle('rank position comparison');
            status_tab = zeros(probe_set_num, tot_query_times, tot_query_times);
            for qti=0:tot_query_times-1
                rank_x = cell2mat(rank_detail(:,max(2*qti,1)));

                for qtj = qti+1:tot_query_times
                    rank_y = cell2mat(rank_detail(:,2*qtj));

                    status_tab(:,qtj,qti+1) = rank_y - rank_x;
                    subplot(tot_query_times, tot_query_times, qti*tot_query_times+qtj);
                    for probek=1:probe_set_num
                        if status_tab(probek,qtj,qti+1)>0
                            plot_red = plot(rank_x(probek), rank_y(probek), 'o','Color','red'); hold on;
                        elseif status_tab(probek,qtj,qti+1)<0
                            plot_green = plot(rank_x(probek), rank_y(probek), 'o','Color','green'); hold on;
                        else
                            plot_blue = plot(rank_x(probek), rank_y(probek), 'o','Color','blue'); hold on;
                        end

                    end
                    better_num = sum(status_tab(:,qtj,qti+1)<0);
                    worse_num = sum(status_tab(:,qtj,qti+1)>0);
                    same_num = sum(status_tab(:,qtj,qti+1)==0);
                    max_rank = max(max(rank_x), max(rank_y));
                    line([1,max_rank],[1,max_rank],'Color','k','LineStyle','--'); grid on; 
                    axis_range = [1 50:50:max_rank];
                    set(gca,'xtick', axis_range);
                    set(gca,'xticklabel',num2cell(axis_range));
                    set(gca,'ytick', axis_range);
                    set(gca,'yticklabel',num2cell(axis_range));
                    if exist('plot_red','var') && exist('plot_green','var') && exist('plot_blue','var')
                        legend([plot_red(1), plot_green(1), plot_blue(1)], ...
                            {'worse','better','same'},'Location','northwest');
                    elseif exist('plot_red','var') && exist('plot_blue','var')
                        legend([plot_red(1), plot_blue(1)], ...
                            {'worse','same'},'Location','northwest');
                    elseif exist('plot_red','var') && exist('plot_green','var')
                        legend([plot_red(1), plot_green(1)], ...
                            {'worse','better'},'Location','northwest');
                    elseif exist('plot_green','var') && exist('plot_blue','var')
                        legend([plot_red(1), plot_green(1)], ...
                            {'better','same'},'Location','northwest');
                    elseif exist('plot_red','var')
                        legend(plot_red(1), {'worse'},'Location','northwest');
                    elseif exist('plot_green','var')
                        legend(plot_green(1), {'better'},'Location','northwest');
                    elseif exist('plot_blue','var')
                        legend(plot_blue(1), {'same'},'Location','northwest');
                    end

                    if qti==0
                        xlabel('initial rank');
                    else
                        xlabel(sprintf('qt=%d rank',qti));
                    end
                    ylabel(sprintf('qt=%d rank',qtj));
                    title(sprintf('#better:%d(%.0f%%) | #worse:%d(%.0f%%) | #same:%d(%.0f%%)', ...
                        better_num, 100*better_num/probe_set_num, ...
                        worse_num, 100*worse_num/probe_set_num, ...
                        same_num, 100*same_num/probe_set_num));
                end
            end
            saveas(hfig, [result_mat_file(1:end-4), '-fb_details_pairwise.fig']);
            if ~show_figure_flag, close(hfig); end

            %% feedback details - figure version
            hfig = figure;
            edges = 0:10:max_rank;
            init_rank = cell2mat(rank_detail(:,1));
            for qt=1:tot_query_times
                subplot(3,tot_query_times,qt);
                rank = init_rank(status_tab(:,qt,1)<0);
                histogram(rank,edges,'FaceColor', 'green');
                xlabel('initial rank'); ylabel('#probes'); title(sprintf('better (qt=%d)',qt));


                subplot(3,tot_query_times,qt+tot_query_times);
                rank = init_rank(status_tab(:,qt,1)>0);
                histogram(rank,edges,'FaceColor', 'red');
                xlabel('initial rank'); ylabel('#probes'); title(sprintf('worse (qt=%d)',qt));


                subplot(3,tot_query_times,qt+2*tot_query_times);
                rank = init_rank(status_tab(:,qt,1)==0);
                histogram(rank,edges,'FaceColor', 'blue');
                xlabel('initial rank'); ylabel('#probes'); title(sprintf('same (qt=%d)',qt));
            end
            saveas(hfig, [result_mat_file(1:end-4), '-fb_details_histogram.fig']);
            if ~show_figure_flag, close(hfig); end
        end
        
end