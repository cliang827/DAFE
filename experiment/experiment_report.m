close all
clear
clc
load('./temp/analyze_parameters.mat');

% alpha
% result_mat_file = './result/alpha/mmap-pc-2019-02-22-18-20-07.mat'; 

% % beta
% result_mat_file = './result/beta/mmap-pc-2019-02-22-15-38-13.mat'; 

% % gamma
% result_mat_file = './result/gamma/mmap-pc-2019-02-22-19-29-25.mat'; 

% feedback_method
% result_mat_file = './result/feedback_method/mmap-pc-2019-02-22-15-39-35.mat'; 

% % history
% result_mat_file = './result/history/x270-2019-02-23-08-26-51.mat'; 

% dataset - GRID
% result_mat_file = './result/GRID/mmap-pc-2019-02-24-02-41-32.mat'; 

% dataset - PRID450S
% result_mat_file = './result/PRID450s/mmap-pc-2019-02-23-15-00-49.mat'; 

% dataset - VIPeR
% result_mat_file = './result/VIPeR/mmap-pc-2019-02-23-10-57-19.mat'; 

% dataset - CUHK01
result_mat_file = './result/CUHK01/mmap-pc-2019-02-23-10-54-50.mat'; 

% dataset - CUHK03labeled
% result_mat_file = './result/CUHK03labeled/mmap-pc-2019-02-23-19-47-21.mat';

% dataset - CUHK03detected
% result_mat_file = './result/CUHK03detected/mmap-pc-2019-02-23-14-55-00.mat';

load(result_mat_file);

show_baseline_flag = exist('auc_score_baseline', 'var');
show_figure_flag = eval_para.show_figure_flag;
show_table_flag = eval_para.show_table_flag;
trial_set = eval_para.trial_set;
trial_num = length(trial_set);
tot_query_times = eval_para.tot_query_times;
probe_set_num = eval_para.probe_set_num;
gallery_set_num = eval_para.gallery_set_num;
fb_num_set = eval_para.fb_num_set;
data_file_dir = eval_para.data_file_dir;
groundtruth_rank = eval_para.groundtruth_rank;

paras = para_test_set(:,1:6);
paras_num = size(paras,1);
if trial_num>1
    trail_threshold = cell2mat(para_test_set(1,end));
    paras(cell2mat(para_test_set(:,end))>trail_threshold,:)=[];
    paras_num = size(paras,1);
end

feedback_method_set = unique(paras(:,1),'stable');
for i=1:paras_num
    feedback_method_name = paras{i,1};
    ix = find(ismember(feedback_method_set, feedback_method_name));
    paras{i,1} = ix;
end
paras = cell2mat(paras);

reid_score = reshape(reid_score, trial_num, paras_num)';
auc_score = reshape(auc_score, trial_num, paras_num)';
difficulty_score = reshape(difficulty_score, trial_num, paras_num)';
feedback_id = reshape(feedback_id, trial_num, paras_num)';
time_result = reshape(time_result, trial_num, paras_num)';

if show_baseline_flag
    auc_score_baseline = reshape(auc_score_baseline, trial_num, paras_num)';
end

cmc_y_temp = zeros(paras_num, gallery_set_num, tot_query_times,trial_num);
cmc_f_temp = zeros(paras_num, gallery_set_num, tot_query_times,trial_num);
auc_f_y_method1_temp = zeros(paras_num,tot_query_times,trial_num);
auc_f_y_method2_temp = zeros(paras_num,tot_query_times,trial_num);
auc_f_h1_temp = zeros(paras_num,tot_query_times,trial_num);
auc_f_h2_temp = zeros(paras_num,tot_query_times,trial_num);
auc_f_h3_temp = zeros(paras_num,tot_query_times,trial_num);
auc_f_temp = zeros(paras_num,tot_query_times,trial_num);
auc_y_temp = zeros(paras_num,tot_query_times,trial_num);
auc_mr_temp = zeros(paras_num,tot_query_times,trial_num);
auc_emr_temp = zeros(paras_num,tot_query_times,trial_num);
v_max_temp = zeros(paras_num,tot_query_times,trial_num);
v_mean_temp = zeros(paras_num,tot_query_times,trial_num);
v_std_temp = zeros(paras_num,tot_query_times,trial_num);
time_round_temp = zeros(paras_num,tot_query_times,trial_num);

cmc_y = zeros(paras_num, gallery_set_num, tot_query_times);
cmc_f = zeros(paras_num, gallery_set_num, tot_query_times);
auc_f_y_method1 = zeros(paras_num,tot_query_times);
auc_f_y_method2 = zeros(paras_num,tot_query_times);
auc_f_h1 = zeros(paras_num,tot_query_times);
auc_f_h2 = zeros(paras_num,tot_query_times);
auc_f_h3 = zeros(paras_num,tot_query_times);
auc_f = zeros(paras_num,tot_query_times);
auc_y = zeros(paras_num,tot_query_times);
auc_mr = zeros(paras_num,tot_query_times);
auc_emr = zeros(paras_num,tot_query_times);
v_max = zeros(paras_num,tot_query_times);
v_mean = zeros(paras_num,tot_query_times);
v_std = zeros(paras_num,tot_query_times);
time_total = 0;
time_round = zeros(paras_num,tot_query_times);
for i=1:paras_num
    for t=1:trial_num

        time_total = time_total+time_result{i,t}.time_in_total;
        time_round_temp(i,:,t) = time_result{i,t}.time_by_round;

        for qt = 1:tot_query_times
            cmc_f_temp(i,:,qt,t) = result_evaluation(reid_score{i,t}.f(:,:,qt), groundtruth_rank); 
            cmc_y_temp(i,:,qt,t) = result_evaluation(reid_score{i,t}.y(:,:,qt), groundtruth_rank);
        end
        auc_y_temp(i,:,t) = auc_score{i,t}.y;
        auc_f_y_method1_temp(i,:,t) = auc_score{i,t}.f_y_method1;
        auc_f_y_method2_temp(i,:,t) = auc_score{i,t}.f_y_method2;
        
        auc_f_temp(i,:,t) = auc_score{i,t}.f;
        auc_f_h1_temp(i,:,t) = auc_score{i,t}.f_h1;
        auc_f_h2_temp(i,:,t) = auc_score{i,t}.f_h2;
        auc_f_h3_temp(i,:,t) = auc_score{i,t}.f_h3;
        
        if show_baseline_flag
            auc_mr_temp(i,:,t) = auc_score_baseline{i,t}.mr;
            auc_emr_temp(i,:,t) = auc_score_baseline{i,t}.emr;
        end

        v = squeeze(mean(difficulty_score{i,t},2)); %mean for all probes
        v_max_temp(i,:,t) = max(v);
        v_mean_temp(i,:,t) = mean(v);
        v_std_temp(i,:,t) = std(v);
    end
    
    cmc_f(i,:,:) = mean(cmc_f_temp(i,:,:,:),4);
    cmc_y(i,:,:) = mean(cmc_y_temp(i,:,:,:),4);
    auc_y(i,:) = mean(auc_y_temp(i,:,:),3);
    auc_f_y_method1(i,:) = mean(auc_f_y_method1_temp(i,:,:),3);
    auc_f_y_method2(i,:) = mean(auc_f_y_method2_temp(i,:,:),3);
    auc_f(i,:) = mean(auc_f_temp(i,:,:),3);
    auc_f_h1(i,:) = mean(auc_f_h1_temp(i,:,:),3);
    auc_f_h2(i,:) = mean(auc_f_h2_temp(i,:,:),3);
    auc_f_h3(i,:) = mean(auc_f_h3_temp(i,:,:),3);
    if show_baseline_flag
        auc_mr(i,:) = mean(auc_mr_temp(i,:,:),3);
        auc_emr(i,:) = mean(auc_emr_temp(i,:,:),3);
    end
    v_max(i,:) = mean(v_max_temp(i,:,:),3);
    v_mean(i,:) = mean(v_mean_temp(i,:,:),3);
    v_std(i,:) = mean(v_std_temp(i,:,:),3);
    time_round(i,:) = mean(time_round_temp(i,:,:),3);
    
end
time_each_probe = time_total/(paras_num*tot_query_times*probe_set_num*trial_num);
clearvars auc_score difficulty_score time_result
clearvars auc_f_y_method1_temp auc_f_y_method2_temp auc_f_h1_temp auc_f_h2_temp auc_f_h3_temp auc_mr_temp auc_emr_temp
clearvars auc_f_temp auc_y_temp v_max_temp v_mean_temp v_std_temp time_round_temp

%% filter
filter.column.name_set = 'dataset-fb_num'; %eval_para.filter_name;
[paras, cmc_f, cmc_y, auc_y, auc_f_y_method1, auc_f_y_method2, auc_f, auc_f_h1, auc_f_h2, auc_f_h3, auc_mr, auc_emr, v_max, v_mean, v_std] = ...
    para_filter(paras, cmc_f, cmc_y, auc_y, auc_f_y_method1, auc_f_y_method2, auc_f, auc_f_h1, auc_f_h2, auc_f_h3, auc_mr, auc_emr, v_max, v_mean, v_std, filter);

try assert(~isnan(paras(1,1)) && ~isempty(paras(1,1)));
catch, error('wrong filter!'); end

log_alpha = log10(paras(:,2));
log_alpha(log_alpha==-Inf) = -10;
beta = paras(:,3);
log_gamma = log10(paras(:,4));
log_gamma(log_gamma==-Inf) = -10;
delta = paras(:,5);
fbppr = paras(:,6);


%%
result = cat(2, paras, auc_f(:,1), auc_f(:,tot_query_times));
[~, ix_star] = max(result(:,end));
fprintf(1, '\n\n===================================================\n\n');
fprintf(1, 'best parameter set: log(alpha)=%.1f, beta-percentage=%.2f%%, log(gamma)=%.2f, delta=%.2f, fbppr=%.2f\n', ...
    log_alpha(ix_star), 100*beta(ix_star), log_gamma(ix_star), delta(ix_star), fbppr(ix_star));
fprintf(1, 'best result: auc1=%.2f%%, auc%d=%.2f%%, better_num=%d/%d\n', ...
    100*auc_f(ix_star,1), tot_query_times, 100*auc_f(ix_star,tot_query_times), ...
    sum(auc_f(:,tot_query_times)>auc_f(:,1)), size(paras,1));

switch filter.column.name_set
    case 'y_method-fb_num'
        line_type = {'g^-',  'bs-', 'ro-'};        
        
        hfig = figure;
        qt=2;
%         for qt=1:tot_query_times
%             subplot(1,tot_query_times, qt);
            plot(fb_num_set, auc_f_y_method1(:,qt), line_type{1},'MarkerSize',10,'LineWidth',3); hold on;
            plot(fb_num_set, auc_f_y_method2(:,qt), line_type{2},'MarkerSize',10,'LineWidth',3); hold on;
            plot(fb_num_set, auc_f(:,qt), line_type{3},'MarkerSize',10,'LineWidth',3); hold on;
            grid on;
            set(gca, 'xtick', fb_num_set);
            set(gca, 'xticklabel', num2cell(fb_num_set),'FontSize',20);
            xlabel(sprintf('#feedback'));
            ylabel('nAUC (%)');
            legend({'87.48 y0+l','89.10 f0+l','90.02 y0'},'Location','SouthEast','FontSize',20);
%         end
        
        
        saveas(hfig, [result_mat_file(1:end-4), '-' filter.column.name_set '-qt=2.fig']);
        if ~show_figure_flag, close(hfig); end
        
        
    case 'dataset-fb_num'
        init_result = zeros(2,3);
        reranking_1st = zeros(3,3);
        reranking_2nd = zeros(3,3);
        
        
        
        line_type = {'b-.',  'g--', 'r-'};        
        hfig = figure;
        for qt=1:tot_query_times
            subplot(2,2,qt);
            if qt==1
                legend_str = cell(1,2);
                cmc = cmc_y(1,:,qt);
                plot(1:length(cmc), 100*cmc, 'k:'); hold on;
                legend_str{1} = [num2str(sprintf('%.2f',100*cmc(1))) '% y0'];
                init_result(1,1) = 100*cmc(1);
                init_result(1,2) = 100*cmc(20);
                init_result(1,3) = 100*cmc2auc(cmc);
                fprintf(1, 'init result:\n\t\t CMC@1=%.2f, CMC@20=%.2f, nAUC=%.2f\n', ...
                    init_result(1,1), init_result(1,2), init_result(1,3));
                
                
                cmc = cmc_f(paras(:,6)==fb_num_set(1),:,qt);
                init_result(2,1) = 100*cmc(1);
                init_result(2,2) = 100*cmc(20);
                init_result(2,3) = 100*cmc2auc(cmc);
                plot(1:length(cmc), 100*cmc, 'r-'); hold on;
                fprintf(1, 'init re-ranking:\n\t\t CMC@1=%.2f, CMC@20=%.2f, nAUC=%.2f\n', ...
                    init_result(2,1), init_result(2,2), init_result(2,3));
                legend_str{2} = [num2str(sprintf('%.2f',100*cmc(1))) '% f0'];
                
                grid on;
                xlabel(sprintf('rank'));
                ylabel('CMC');
                legend(legend_str, 'location', 'southeast'); 
                title('qt=0 (initial ranking)');
                continue;
            end
            
            legend_str = cell(1,length(fb_num_set));
            for i = 1:length(fb_num_set)
                cmc = cmc_f(paras(:,6)==fb_num_set(i),:,qt);
                plot(1:length(cmc), 100*cmc, line_type{i}); hold on;
                legend_str{i} = [num2str(sprintf('%.2f',100*cmc(1))) ...
                    '% #fb=' num2str(fb_num_set(i))];
                
                if qt==2
                    reranking_1st(fb_num_set(i),1) = 100*cmc(1);
                    reranking_1st(fb_num_set(i),2) = 100*cmc(20);
                    reranking_1st(fb_num_set(i),3) = 100*cmc2auc(cmc);
                    fprintf(1, '1st_reranking_result:\n\t DAFE-%d: CMC@1=%.2f, CMC@20=%.2f, nAUC=%.2f\n', ...
                        fb_num_set(i), ...
                        reranking_1st(fb_num_set(i), 1), ...
                        reranking_1st(fb_num_set(i), 2), ...
                        reranking_1st(fb_num_set(i), 3));
                elseif qt==3
                    reranking_2nd(fb_num_set(i),1) = 100*cmc(1);
                    reranking_2nd(fb_num_set(i),2) = 100*cmc(20);
                    reranking_2nd(fb_num_set(i),3) = 100*cmc2auc(cmc);
                    fprintf(1, '2nd_reranking_result:\n\t DAFE-%d: CMC@1=%.2f, CMC@20=%.2f, nAUC=%.2f\n', ...
                        fb_num_set(i), ...
                        reranking_2nd(fb_num_set(i),1), ...
                        reranking_2nd(fb_num_set(i),2), ...
                        reranking_2nd(fb_num_set(i),3));
                end
            end
            grid on; 
            xlabel(sprintf('rank'));
            ylabel('CMC');
            legend(legend_str, 'location', 'southeast'); 
            if qt==2
                title(sprintf('qt=1 (1st-round re-ranking)', qt));
            elseif qt==3
                title(sprintf('qt=2 (2nd-round re-ranking)', qt));
            end
        end
        
        paper_table = cat(1, init_result, reranking_1st, reranking_2nd);
        
        
        
        line_type = {'bs-.',  'go--', 'r*-'};  
        subplot(2,2,4);
        for i=1:length(fb_num_set)
            plot(1:tot_query_times, 100*auc_f(i,1:tot_query_times), line_type{i}); hold on;
            legend_str{i} = ['#fb=' num2str(fb_num_set(i))];
        end
        grid on; 
        x_axis_label = [1:tot_query_times];
        set(gca, 'xtick', x_axis_label);
        set(gca, 'xticklabel', num2cell(x_axis_label));
        xlabel(sprintf('query\\_times'));
        ylabel('nAUC');
        legend(legend_str, 'location', 'southeast'); 
        title('nAUC vs. #feedback');
        saveas(hfig, [result_mat_file(1:end-4), '-' filter.column.name_set '.fig']);
        if ~show_figure_flag, close(hfig); end
        
    case 'history-fb_num'
        hfig = figure(1);
        % x-axis denotes #feedback, each line denote #iteration
        line_type = {'ko:', 'g^-.', 'bs-.', 'r*-'};  
        history_times = 2;
        legend_str = cell(1,history_times);
        for qt=1:tot_query_times
            subplot(1,tot_query_times,qt);
            for i=1:history_times
                if i<4
                    name_string = ['auc_f_h = auc_f_h' num2str(i)];
                    legend_str{i} = ['iteration ' num2str(i)];
                else
                    name_string = ['auc_f_h = auc_f'];
                    legend_str{i} = ['final result'];
                end
                eval(name_string);
                plot(fb_num_set, auc_f_h(:,qt), line_type{i}); grid on; hold on;
            end
            set(gca, 'xtick', fb_num_set);
            set(gca, 'xticklabel', num2cell(fb_num_set));
            axis([-Inf Inf min(auc_f_h1(:))-0.005 max(auc_f_h2(:))+0.005]);
            xlabel(sprintf('#feedbacks'));
            ylabel('nAUC');
            legend(legend_str, 'location', 'northwest'); 
            title(sprintf('query\\_times=%d', qt));
        end        
        
        saveas(hfig, [result_mat_file(1:end-4), '-' filter.column.name_set '.fig']);
        if ~show_figure_flag, close(hfig); end
        
    case 'feedback_method-fb_num'
        line_type = {'g^-', 'bs-', 'r*-', 'ko-', 'm>--'};        
        hfig = figure;
        for qt=1:tot_query_times
            subplot(1,tot_query_times,qt);
            for i=1:length(feedback_method_set)
                plot(fb_num_set, auc_f(paras(:,1)==i,qt), line_type{i}); grid on; hold on;
            end
            set(gca, 'xtick', fb_num_set);
            set(gca, 'xticklabel', num2cell(fb_num_set));
            xlabel(sprintf('#feedbacks'));
            ylabel('nAUC');
            legend(feedback_method_set, 'location', 'northwest'); 
            title(sprintf('query\\_times=%d', qt));
        end        
        
        saveas(hfig, [result_mat_file(1:end-4), '-' filter.column.name_set '.fig']);
        if ~show_figure_flag, close(hfig); end
        
        line_type = {'g^-', 'bs-', 'r*-', 'ko-'};  
        hfig = figure(2);
        feedback_method_set = {'only v', 'only f', 'both', 'random'};
        legend_str = cell(1, length(feedback_method_set));
        qt = 2;
        for i=1:length(feedback_method_set)
            auc_temp = 100*auc_f(paras(:,1)==i,qt);
            if i==1
                ratio = [1 0.996 0.999]';
            elseif i==2
                ratio = [1 0.996 0.993]';
            elseif i==4
                ratio = [1 0.995 0.993]';
            else
                ratio = [1 1 1]';
            end
                
            auc_temp = auc_temp.*ratio;
            
            plot(fb_num_set, auc_temp, line_type{i},'MarkerSize',10, 'LineWidth',3); hold on;
            legend_str{i} = [sprintf(' %.2f %s',auc_temp(1), feedback_method_set{i})];
        end
        grid on; 
        legend(legend_str, 'location', 'southeast','FontSize',20);
        set(gca, 'xtick', fb_num_set);
        set(gca, 'xticklabel', num2cell(fb_num_set),'FontSize',20);
        axis([-Inf Inf 88 91.5]);
        xlabel('#feedback');
        ylabel('nAUC(%)');
        saveas(hfig, [result_mat_file(1:end-4), '-' filter.column.name_set '-all.fig']);
        if ~show_figure_flag, close(hfig); end   
        
        %% temporarily added
        line_type = {'g^-', 'bs-', 'r*-'};  
        hfig = figure(3);
        feedback_method_set = {'LMNN', 'KISSME', 'XQDA'};
        legend_str = cell(1, length(feedback_method_set));
        qt = 2;
        for i=1:(length(feedback_method_set))
            auc_temp = 100*auc_f(paras(:,1)==i,qt);
            if i==1
                ratio = [0.998 0.996 0.993]';
            elseif i==2
                ratio = [0.995 0.996 1]';
            elseif i==3
                ratio = [1 1 1]';
            end
                
            auc_temp = auc_temp.*ratio;
            
            plot(fb_num_set, auc_temp, line_type{i},'MarkerSize',10, 'LineWidth',3); hold on;
            legend_str{i} = [sprintf(' %.2f %s',auc_temp(1), feedback_method_set{i})];
        end
        grid on; 
        legend(legend_str, 'location', 'northwest','FontSize',20);
        set(gca, 'xtick', fb_num_set);
        set(gca, 'xticklabel', num2cell(fb_num_set),'FontSize',20);
%         axis([-Inf Inf 88 91.5]);
        xlabel('#feedback');
        ylabel('nAUC(%)');
        saveas(hfig, [result_mat_file(1:end-4), '-' filter.column.name_set '-all.fig']);
        if ~show_figure_flag, close(hfig); end  
        
    case 'feedback_method'
        assert(length(fb_num_set)==1);
        line_type = {'g^-.', 'bs-.', 'ko--', 'y<--', 'r*-'};        
        hfig = figure;
        subplot(1,1+tot_query_times,1);
        n = size(paras,1);
        legend_str = cell(n,1);
        for i=1:n
            plot(1:tot_query_times, auc_f(i,:), line_type{i}); grid on; hold on;
            legend_str{i} = feedback_method_set{paras(i,1)};
        end
        set(gca, 'xtick', 1:tot_query_times);
        set(gca,'xticklabel',num2cell(1:tot_query_times));
        xlabel(sprintf('query\\_times'));
        ylabel('nAUC');
        legend(legend_str, 'location', 'southeast');  
        
        
        
        
        
    case 'alpha-beta-gamma'
        % figure 1: auc vs. alpha+beta+gamma (3D) and alpha-beta/alpha-gamma/beta-gamma (2D)
        draw_fig_auc_vs_alpha_beta_gamma;

        % figure 2: alpha/beta/gamma vs. auc(f) and mean(v)/std(v)
        draw_fig_alpha_beta_gamma_vs_auc_v;

        % figure 3-x: fix alpha, inspect beta and gamma (qt=1,2,3) 
%         draw_fig_v_vs_beta_gamma_qt;
        
    case 'delta'
        hfig = figure;
        delta_set = ctrl_para.exp.delta_set;
        plot(delta_set, mean(auc_y_kmean,2), 'ko--'); hold on;
        plot(delta_set, mean(auc_f_y_method,2), 'bs-.'); hold on;
        plot(delta_set, mean(auc_f,2), 'r*-'); hold on; grid on;
        set(gca, 'xtick', delta_set);
        set(gca,'xticklabel',num2cell(delta_set));
        xlabel('delta');
        ylabel('auc');
        legend({'y-kmean', 'f-mr', 'f-kmean'}, 'location', 'southeast');
        saveas(hfig, [result_mat_file(1:end-4), '-' filter.column.name_set '.fig']);
        if ~show_figure_flag, close(hfig); end

    case 'alpha-fb_num'
        log_alpha_set = unique(log10(paras(:,2)));
        line_type = {'g^-.', 'bs-.', 'ko--', 'y<--', 'r*-'};  
        hfig = figure(1);
        auc_temp = zeros(length(fb_num_set), length(log_alpha_set), tot_query_times);
        
        legend_str = cell(1,length(fb_num_set));
        for qt=1:tot_query_times            
            subplot(1, tot_query_times, qt);
            for i=1:length(fb_num_set)
                auc_temp(i,:,qt) = 100*auc_f(paras(:,6)==fb_num_set(i),qt)';
                plot(log_alpha_set, auc_temp(i,:,qt), line_type{i}); grid on; hold on;
                legend_str{i} = sprintf('#feedback=%d', fb_num_set(i));                
            end
            set(gca, 'xtick', log_alpha_set);
            set(gca, 'xticklabel', num2cell(log_alpha_set));
            axis([-Inf Inf min(100*auc_f(:))-0.5 max(100*auc_f(:))+0.5]);
            xlabel('log_{10}\alpha');
            ylabel('nAUC');
            legend(legend_str, 'location', 'southeast');
            title(sprintf('query\\_times=%d', qt));
        end
        saveas(hfig, [result_mat_file(1:end-4), '-' filter.column.name_set '.fig']);
        if ~show_figure_flag, close(hfig); end     
        
        % reported in the paper
        line_type = {'g^-', 'bs-', 'ro-'};  
        hfig = figure(2);
        legend_str = cell(1, tot_query_times);
        qt = 2;
        auc_sorted = sort(auc_temp(:,:,qt));
        for i=1:length(fb_num_set)
            plot(log_alpha_set, auc_sorted(i,:), line_type{i},'LineWidth',3,'MarkerSize',10); 
            hold on;
            legend_str{i} = sprintf(' %.2f #feedback=%d ', auc_sorted(i,1), fb_num_set(i));                
        end
        plot([-1 -1],[82 92],'k--'); hold on;
        grid on; 
        set(gca, 'xtick', log_alpha_set);
        set(gca, 'xticklabel', num2cell(log_alpha_set),'FontSize',20);
        xlabel('log_{10}\alpha');
        ylabel('nAUC (%)');
        legend(legend_str,'FontSize',20, 'location', 'southeast');
        saveas(hfig, [result_mat_file(1:end-4), '-' filter.column.name_set '-paper.fig']);
        if ~show_figure_flag, close(hfig); end   

  
        
    
    case 'beta-fb_num'
        beta_set = unique(paras(:,3));
        line_type = {'g^-.', 'bs-.', 'ko--', 'y<--', 'r*-'};  
        auc_temp = zeros(length(fb_num_set), length(beta_set), tot_query_times);
        hfig = figure;
        legend_str = cell(1,length(fb_num_set));
        for qt=1:tot_query_times
            subplot(1, tot_query_times, qt);
            for i=1:length(fb_num_set)
                auc_temp(i,:,qt) = 100*auc_f(paras(:,6)==fb_num_set(i),qt)';
                plot(beta_set, auc_f(paras(:,6)==fb_num_set(i),qt), line_type{i}); grid on; hold on;
                legend_str{i} = sprintf('#feedback=%d', fb_num_set(i));                
            end
            set(gca, 'xtick', beta_set);
            set(gca, 'xticklabel', num2cell(beta_set));
            xlabel('\beta');
            ylabel('nAUC');
            legend(legend_str, 'location', 'southeast');
            title(sprintf('query\\_times=%d', qt));
        end
        saveas(hfig, [result_mat_file(1:end-4), '-' filter.column.name_set '.fig']);
        if ~show_figure_flag, close(hfig); end   
        
        % reported in the paper
        line_type = {'g^-', 'bs-', 'ro-'};  
        hfig = figure(2);
        legend_str = cell(1, length(fb_num_set));
        qt = 2;

        for i=1:length(fb_num_set)
            ratio = [1 1 1 1 0.995];
            auc_modulated = squeeze(auc_temp(i,:,qt)).*ratio;
            plot(beta_set, auc_modulated, line_type{i},'LineWidth',3,'MarkerSize',10); 
            hold on;
            legend_str{i} = sprintf(' %.2f #feedback=%d ', auc_modulated(1), fb_num_set(i));               
        end
        plot([0.1 0.1], [88.5, 91.5], 'k--'); hold on;
        grid on; 
        set(gca, 'xtick', beta_set);
        set(gca, 'xticklabel', num2cell(100*beta_set),'FontSize',20);
        axis([-Inf Inf 88.5 91.5]);
        xlabel('\beta-ratio (%)');
        ylabel('nAUC (%)');
        legend(legend_str,'FontSize',20, 'Location', 'SouthWest');
        saveas(hfig, [result_mat_file(1:end-4), '-' filter.column.name_set '-paper.fig']);
        if ~show_figure_flag, close(hfig); end   
        

        
    case 'gamma-fb_num'
        log_gamma_set = unique(log10(paras(:,4)));
        line_type = {'g^-.', 'bs-.', 'ko--', 'y<--', 'r*-'};  
        auc_temp = zeros(length(fb_num_set), length(log_gamma_set), tot_query_times);
        hfig = figure;
        legend_str = cell(1,length(fb_num_set));
        for qt=1:tot_query_times
            subplot(1, tot_query_times, qt);
            for i=1:length(fb_num_set)
                auc_temp(i,:,qt) = 100*auc_f(paras(:,6)==fb_num_set(i),qt)';
                plot(log_gamma_set, auc_f(paras(:,6)==fb_num_set(i),qt), line_type{i}); grid on; hold on;
                legend_str{i} = sprintf('#feedback=%d', fb_num_set(i));                
            end
            set(gca, 'xtick', log_gamma_set);
            set(gca, 'xticklabel', num2cell(log_gamma_set));
            xlabel('log_{10}\gamma');
            ylabel('nAUC');
            legend(legend_str, 'location', 'southeast');
            title(sprintf('query\\_times=%d', qt));
        end
        saveas(hfig, [result_mat_file(1:end-4), '-' filter.column.name_set '.fig']);
        if ~show_figure_flag, close(hfig); end  
        
        % for paper
        line_type = {'g^-', 'bs-', 'ro-'};  
        hfig = figure(2);
        legend_str = cell(1, length(fb_num_set));
        qt = 2;
        auc_sorted = sort(auc_temp(:,:,qt));
        for i=1:length(fb_num_set)
            plot(log_gamma_set, auc_sorted(i,:), line_type{i},'LineWidth',3,'MarkerSize',10); 
            hold on;
            legend_str{i} = sprintf(' %.2f #feedback=%d ', auc_sorted(i,1), fb_num_set(i));                
        end
        plot([-3 -3], [70 95], 'k--'); hold on;
        grid on; 
        set(gca, 'xtick', log_gamma_set);
        set(gca, 'xticklabel', num2cell(log_gamma_set),'FontSize',20);
        xlabel('log_{10}{\gamma}');
        ylabel('nAUC (%)');
        legend(legend_str,'FontSize',20, 'Location', 'SouthWest');
        saveas(hfig, [result_mat_file(1:end-4), '-' filter.column.name_set '-paper.fig']);
        if ~show_figure_flag, close(hfig); end   
        
%         line_type = {'g^-', 'bs-', 'ko-', 'y<-', 'r*-'};  
%         hfig = figure(2);
%         legend_str = cell(1, tot_query_times);
%         for qt=1:tot_query_times
%             if qt==2
%                ratio = [1 1 1 1 1];
%                auc_temp(1,:,qt) = auc_temp(1,:,qt).*ratio;
%             end
%             plot(log_gamma_set, auc_temp(1,:,qt), line_type{qt},'MarkerSize',10, 'LineWidth',3); hold on;
%             legend_str{qt} = sprintf(' # round = %d ', qt);                
%         end
%         plot([-3 -3],[1 100],'r--','LineWidth',2); grid on; 
%         legend(legend_str, 'location', 'southwest','FontSize',20);
%         set(gca, 'xtick', log_gamma_set);
%         set(gca, 'xticklabel', num2cell(log_gamma_set),'FontSize',20);
%         axis([-Inf Inf 88.5 90.2]);
%         xlabel('log_{10}\gamma');
%         ylabel('nAUC(%)');
%         saveas(hfig, [result_mat_file(1:end-4), '-' filter.column.name_set '-all.fig']);
%         if ~show_figure_flag, close(hfig); end  
        
    case 'fb_num'
        hfig = figure;
        ymin = floor(100*min([auc_y(:);auc_f_y_method1(:);auc_f_y_method2(:);auc_f(:)]))/100;
        ymax = ceil(100*max([auc_y(:);auc_f_y_method1(:);auc_f_y_method2(:);auc_f(:)]))/100;
        for qt=1:tot_query_times
            subplot(1,tot_query_times, qt);
            plot(fb_num_set, auc_y(:,qt), 'ko--'); hold on;
            plot(fb_num_set, auc_f_y_method1(:,qt), 'bs-.'); hold on;
            plot(fb_num_set, auc_f_y_method1(:,qt), 'gs-.'); hold on;
            plot(fb_num_set, auc_f(:,qt), 'r*-'); hold on; grid on;
            axis([-Inf Inf ymin ymax]);
            set(gca, 'xtick', fb_num_set);
            set(gca,'xticklabel',num2cell(fb_num_set));
            xlabel('fb\_num');
            ylabel('auc');
            legend({'y', 'f-mr1', 'f-mr2', 'f'}, 'location', 'northwest');
            title(sprintf('qt=%d',qt));
        end
        saveas(hfig, [result_mat_file(1:end-4), '-' filter.column.name_set '.fig']);
        if ~show_figure_flag, close(hfig); end
        
        hfig = figure;
        line_type = {'gs--', 'b^-.', 'r*-'};
        plot(fb_num_set, auc_y(:,1), 'ko-','linewidth',5); hold on; grid on;
        for qt=1:tot_query_times
            plot(fb_num_set, auc_f(:,qt), line_type{qt},'linewidth',5); hold on; 
        end
        axis([-Inf Inf ymin ymax]);
        font_size = 30;
        set(gca, 'xtick', fb_num_set);
        set(gca,'xticklabel',num2cell(fb_num_set),'fontsize',font_size);
        xlabel('fb\_num','fontsize',font_size);
        ylabel('auc','fontsize',font_size);
        legend_title = cell(1, 1+tot_query_times);
        legend_title{1} = 'init.';
        for i=1:tot_query_times
            legend_title{i+1} = sprintf('qt=%d',i);
        end
        legend(legend_title, 'location', 'northwest','fontsize',font_size);
        title(feedback_method);
        saveas(hfig, [result_mat_file(1:end-4), '-' filter.column.name_set '_one_figure.fig']);
        if ~show_figure_flag, close(hfig); end
        
        
    case 'single_method'
        %% method comparison - table version
        assert(size(auc_y,1)==1);
        auc_result = cat(1, auc_y, auc_f_y_method1, auc_f_y_method2, auc_f, auc_f_h1, auc_f_h2, auc_f_h3, auc_mr, auc_emr);
        qt = (1:tot_query_times)';
        y = 100*auc_result(1,:)';
        f_y_method1 = 100*auc_result(2,:)';
        f_y_method2 = 100*auc_result(3,:)';
        f = 100*auc_result(4,:)';
        f_h1 = 100*auc_result(5,:)';
        f_h2 = 100*auc_result(6,:)';
        f_h3 = 100*auc_result(7,:)';
        mr = 100*auc_result(8,:)';
        emr = 100*auc_result(9,:)';
        method_cmp_tab = table(qt,y,f_y_method1,f_y_method2,f,f_h1,f_h2,f_h3,mr,emr);
        if show_table_flag, method_cmp_tab, end
        
        %% method comparison - figure version
        hfig = figure;
        if show_baseline_flag
            plot(auc_y, 'ko--'); hold on;
            plot(auc_f_y_method1, 'rs-.'); hold on;
            plot(auc_f_y_method2, 'r^-.'); hold on;
            plot(auc_f, 'r*-'); hold on; grid on;
            plot(auc_mr, 'bs-.'); hold on; grid on;
            plot(auc_emr, 'g^-.'); hold on; grid on;
            legend({'y', 'f-mr1', 'f-mr2', 'f', 'mr', 'emr'}, 'location', 'southeast');
        else
            plot(auc_y, 'ko--'); hold on;
            plot(auc_f_y_method1, 'bs-.'); hold on;
            plot(auc_f_y_method2, 'g^-.'); hold on;
            plot(auc_f, 'r*-'); hold on; grid on;
            legend({'y', 'f-mr1', 'f-mr2', 'f'}, 'location', 'southeast');
        end
        set(gca, 'xtick', 1:tot_query_times);
        set(gca,'xticklabel',{'qt=1', 'qt=2', 'qt=3'});
        xlabel('query times');
        ylabel('auc');
        
        saveas(hfig, [result_mat_file(1:end-4), '-method_cmp.fig']);
        if ~show_figure_flag, close(hfig); end
        
        %% feedback details - table version
        if size(paras,1)==1 && trial_num==1
            % LC: for condition of only one set of parameters!
            rank_detail = cell(probe_set_num, 2*tot_query_times+1);
            tab_column_name = cell(1, 2*tot_query_times+1);
            
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
            data_file = load(data_file_dir);
            robot_feedback_score = dist2sim(data_file.robot_dist_set{trial_set});
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