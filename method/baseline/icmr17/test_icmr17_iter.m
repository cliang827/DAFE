% purpose: test icmr17 on pcm14.mat dataset

clc
clear
close all

pwd_dir = pwd;
if ~strcmp(pwd_dir(end-14:end), 'SalUBM_Feedback')
    error('Hi, please set current working dir as SalUBM_Feedback and Add to Path current dir!');
end
addpath(genpath(fullfile([pwd '/method'])));

dataset = 'cuhk02';
feature = 'pcm14';

data_dir = sprintf('./data/%s/init/%s.mat', dataset, feature);


load(data_dir);
[gallery_num, probe_num] = size(g2p_dist_global);
dist = zeros(gallery_num+1, gallery_num+1);

method_set = {'CPRR','RLRECOM','RKGRAPH','CORGRAPH','RLSIM','RECKNNGRAPH','CONTEXTRR'};
% method_set = {'CONTEXTRR'};
method_num = length(method_set);
test_dir = './method/baseline/icmr17/bin/';
result_dir = './method/baseline/icmr17/result/';
exe_path = [test_dir 'udlf']; 

dist_dir = [test_dir 'dist.txt'];
output_dir = [test_dir 'output.txt'];
config_dir = [test_dir 'test_config.ini'];
new_dist = zeros(gallery_num, probe_num);

for method_index = 1:method_num

    for reranking_times = 1:tot_reranking_times
        tic;
        
        if reranking_times==1
            % for the 1st round reranking (namely 2nd round query)
            % g2p_dist_new is initialized with g2p_dist_global
            g2p_dist_new = g2p_dist_global;
        else
            % from 2rd round reranking (namely 3nd round query) 
            % update g2p_dist_new with last round result
            g2p_dist_new = new_dist;
        end

        for i =1:probe_num
            nchar = fprintf(1, '\t ... testing method: %s | #iter. %d/%d (%.2f%%) | #process. %d/%d (%.2f%%)', ...
                method_set{method_index}, reranking_times, tot_reranking_times, 100*reranking_times/tot_reranking_times, ...
                i, probe_num, 100*i/probe_num);

            % construct dist matrix
            dist(1,2:end) = g2p_dist_new(:,i)';
            dist(2:end,1) = g2p_dist_new(:,i);
            dist(2:end,2:end) = g2g_dist_global;

            % write dist matrix into input.txt
            file = fopen(dist_dir, 'w');
            assert(file~=-1);
            for j = 1:(gallery_num+1)
                for k = 1:(gallery_num+1)
                    fprintf(file, '%f ', dist(j,k));
                end
                fprintf(file,'\n ');
            end
            fclose(file); 

            % execute udlf
            [status, cmdout] = system([exe_path,' ',config_dir,' ',method_set{method_index}]); 
            assert(status==0);

            % collect re-ranking result (in dist format) into result
            result = load(output_dir);
            new_dist(:,i) = result(1,2:end)';

            fprintf(1, repmat('\b', 1, nchar));
        end
        
        reranking_cmc(:, reranking_times) = result_cmc(new_dist, ground_truth, 'ascend');
        reranking_auc(1, reranking_times) = cmc2auc(reranking_cmc(:,reranking_times));
        
        fprintf(1, repmat('\b', 1, nchar));
        fprintf(1, '\t method: %s | #iter. %d/%d (%.2f%%) | run_time: %.2f sec\n', ...
            method_set{method_index}, reranking_times, tot_reranking_times, 100*reranking_times/tot_reranking_times, toc);
    end
    
    save(sprintf('%s%s/%s_%s.mat', result_dir, dataset, feature, method_set{method_index}), ...
        '-v7.3', 'init_cmc', 'init_auc', 'reranking_cmc', 'reranking_auc');
end


% % draw auc figure
% color_tab = rand(method_num,3);
% line_style = {'--',':','-.','-'};
% marker_type = {'o', 's', 'p', 'x', 'h', '^', 'v', '<', '>', '.', '*', '+'};
% 
% 
% fig_handle = figure; % auc figure
% for method_index = 1:method_num
%     load(sprintf('%s%s/%s_%s.mat', result_dir, dataset, feature, method_set{method_index}));
%     plot(1:tot_query_times, [init_auc, reranking_auc], ...
%         'Color', color_tab(method_index,:), ...
%         'LineStyle', line_style{mod(method_index,length(line_style))+1}, ...
%         'Marker', marker_type{mod(method_index,length(marker_type))+1});
%     hold on;
% end
% xlabel('iter.');ylabel('nAUC');
% title('nAUC vs. iter.')
% grid on;
% legend(method_set,'Location','southeast');
% 
% print(fig_handle, sprintf('%s%s/%s_auc.eps', result_dir, dataset, feature), '-depsc2', '-r600');
% savefig(fig_handle, sprintf('%s%s/%s_auc.fig',result_dir, dataset, feature));
% 
% % draw cmc figure
% subplot_array = [...
%     1, 1;   % 1
%     1, 2;   % 2
%     1, 3;   % 3
%     2, 2;   % 4
%     2, 3;   % 5
%     2, 3;   % 6
%     2, 4;	% 7
%     2, 4;   % 8
%     3, 3;   % 9
%     3, 4;   % 10
%     3, 4;   % 11
%     3, 4;   % 12
% ];
% 
% if tot_reranking_times>12
%     error('Hi, too many figures can not be drawn in one figure!');
% end
% 
% subplot_row = subplot_array(tot_reranking_times, 1);
% subplot_col = subplot_array(tot_reranking_times, 2);
% 
% scrsz = get(groot,'ScreenSize');
% fig_handle = figure('Position',[1 scrsz(4) scrsz(3) scrsz(4)]); % full screen
%         
% for reranking_times = 1:tot_reranking_times
%     
%     subplot(subplot_row, subplot_col, reranking_times);
%     result_name = cell(1, 1+method_num);
%     
%     % plot init result
%     plot(1:gallery_num, init_cmc, 'k--'); hold on;
%     result_name{1} = sprintf('init 0: %2.2f%%', 100*init_auc);
% 
%     % plot re-ranking result
%     for method_index = 1:method_num
%         load(sprintf('%s%s/%s_%s.mat', result_dir, dataset, feature, method_set{method_index}));
%         new_cmc = reranking_cmc(:,reranking_times);
%         new_auc = reranking_auc(reranking_times);
% 
%         plot(1:gallery_num, new_cmc, ...
%             'Color', color_tab(method_index,:), ...
%             'LineStyle', line_style{mod(method_index,length(line_style))+1}, ...
%             'Marker', marker_type{mod(method_index,length(marker_type))+1});
%         hold on
% 
%         result_name{1+method_index} = sprintf('iter %d: %2.2f%% - %s', ...
%             reranking_times, 100*new_auc, method_set{method_index});
%     end
% 
%     xlabel('rank');ylabel('CMC');
%     legend(result_name,'Location','southeast');
%     axis([0,gallery_num,0,1]);
%     grid on;
% end
% 
% print(fig_handle, sprintf('%s%s/%s_cmc.eps', result_dir, dataset, feature), '-depsc2', '-r600');
% savefig(fig_handle, sprintf('%s%s/%s_cmc.fig',result_dir, dataset, feature));
