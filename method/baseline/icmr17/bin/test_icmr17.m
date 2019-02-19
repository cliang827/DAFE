% purpose: test icmr17 on pcm14.mat dataset

clc
clear

dataset = 'viper';
feature = 'pcm14';

data_dir = sprintf('./data/%s/init/%s.mat', dataset, feature);
dist_mat_dir = sprintf('./method/baseline/icmr17/result/%s/%s_dist.mat', dataset, feature);
update_data_flag = false;
if ~exist(dist_mat_dir,'file') || update_data_flag
    ctrl_para = load(data_dir);
    init_reid_score = ctrl_para.reid_score;
    g2g_dist_part = ctrl_para.g2g_dist;
    g2p_dist_part = ctrl_para.g2p_dist;
    p2g_dist_part = ctrl_para.p2g_dist;

    [gallery_num,probe_num,part_num] = size(g2p_dist_part);
    g2g_dist = zeros(gallery_num,gallery_num);        % whole body, not body part
    g2p_dist = zeros(gallery_num,probe_num);                  % whole body, not body part
    for k=1:part_num
        g2g_dist = g2g_dist + 0.5*(g2g_dist_part(:,:,k) + g2g_dist_part(:,:,k)');
        g2p_dist = g2p_dist + 0.5*(p2g_dist_part(:,:,k)'+ g2p_dist_part(:,:,k));
    end

    save(dist_mat_dir, '-v7.3', 'g2p_dist', 'g2g_dist', 'init_reid_score');
end

load(dist_mat_dir);
[gallery_num, probe_num] = size(g2p_dist);
dist = zeros(gallery_num+1, gallery_num+1);

% method_set = {'CPRR','RLRECOM','RKGRAPH','CORGRAPH','RLSIM','RECKNNGRAPH','CONTEXTRR'};
method_set = {'CPRR'};
test_dir = './method/baseline/icmr17/bin/';
result_dir = './method/baseline/icmr17/result/';
exe_path = [test_dir 'udlf']; 

dist_dir = [test_dir 'dist.txt'];
output_dir = [test_dir 'output.txt'];
config_dir = [test_dir 'test_config.ini'];

method_num = length(method_set);
for method_index = 1:method_num
    new_dist = zeros(gallery_num, probe_num);
    tic;
    for i =1: probe_num
        nchar = fprintf(1, '\t testing method: %s | #pro. %d/%d (%.2f%%)', ...
            method_set{method_index}, i, probe_num, 100*i/probe_num);

        % construct dist matrix
        dist(1,2:end) = g2p_dist(:,i)';
        dist(2:end,1) = g2p_dist(:,i);
        dist(2:end,2:end) = g2g_dist;
        
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
    fprintf(1, repmat('\b', 1, nchar));
    fprintf(1, '\t method: %s | #process. %d/%d (%.2f%%) | run_time: %.2f\n', ...
        method_set{method_index}, i, probe_num, 100*i/probe_num, toc);

    save(sprintf('%s%s/%s_%s.mat', result_dir, dataset, feature, method_set{method_index}),'new_dist');
end


% draw CMC curve
scrsz = get(groot,'ScreenSize');
figure('Position',[1 scrsz(4) scrsz(3) scrsz(4)]); % full screen
ground_truth = repmat(1:probe_num, gallery_num, 1);
result_name = cell(1, 1+method_num);

% plot init result
init_result = g2p_dist;
init_cmc = result_cmc(init_result, ground_truth, 'ascend');
init_auc = cmc2auc(init_cmc);
result_name{1} = sprintf('init 0: %2.2f%%', 100*init_auc);
plot(1:gallery_num, init_cmc, 'k--'); hold on;

% plot re-ranking result
color_tab = rand(method_num,3);
line_style = {'--',':','-.','-'};
marker_type = {'+', 'o', 's', '*', 'x', '^', 'v', '<', '>', '.', 'p', 'h'};
for method_index = 1:method_num
    load(sprintf('%s%s/%s_%s.mat', result_dir, dataset, feature, method_set{method_index}));
    new_cmc = result_cmc(new_dist, ground_truth, 'ascend');
    new_auc = cmc2auc(new_cmc);
    
    plot(1:gallery_num, new_cmc, ...
        'Color', color_tab(method_index,:), ...
        'LineStyle', line_style{mod(method_index,length(line_style))+1}, ...
        'Marker', marker_type{mod(method_index,length(marker_type))+1});
    hold on
    
    result_name{1+method_index} = sprintf('iter 1: %2.2f%% - %s', 100*new_auc, method_set{method_index});
end

xlabel('rank');ylabel('CMC');
legend(result_name(1:method_index+1),'Location','southeast');
axis([0,gallery_num,0,1]);
set(gca,'fontsize',30);
grid on;

print(sprintf('./method/baseline/icmr17/result/viper/%s_cmc.eps', feature), '-depsc2', '-r600');
print(sprintf('./method/baseline/icmr17/result/viper/%s_cmc.jpg', feature), '-djpeg', '-r666');
