
%% 程序是在linux下调试通过，在Windows下调用exe的函数运行不通过，初步排查是编译器的问题，用vs编译器编译的程序可以调用运行，
%% 而用g++编译器编译的运行不了

clc
clear

method = {'CPRR','RLRECOM','RLSIM','CONTEXTRR','RECKNNGRAPH','RKGRAPH','CORGRAPH'};

exe_path = './udlf'; 
distance_mat_save_path = ['Viper','/','matrices','/','distance.txt'];
reranking_result_path = 'output.txt';
config_file_path = 'config.ini';

load('dis_p2g_all.mat');
load('dis_g2g_all.mat');
reranking_distance = [];
Distance_mat_size = size(dis_p2g_all,2) + 1;
Distance_mat = zeros(Distance_mat_size,Distance_mat_size);
for method_index = 1:length(method)
    for i =1: size(dis_p2g_all,2)
        tic
        %%构建距离矩阵
        Distance_mat(1,1) = 0;
        Distance_mat(1,2:end) = dis_p2g_all(:,i)';
        Distance_mat(2:end,1) = dis_p2g_all(:,i);
        Distance_mat(2:end,2:end) = dis_g2g_all;
        
        %%将距离矩阵写入输入txt矩阵
        file = fopen(distance_mat_save_path,'w');
        assert(file~=-1);
        for ii = 1:size(Distance_mat,1)
            for jj = 1:size(Distance_mat,2)
                fprintf(file,'%f ',Distance_mat(ii,jj));
            end
            fprintf(file,'\n ');
        end
        fclose(file); 
        
        %执行外部算法程序
        [status, cmdout] = system([exe_path,' ',config_file_path,' ',method{method_index}]); 
        t = toc;
        
        if status==0
            fprintf('the UDLF of %s has done! progress %d/%d - %.2f\n',method{method_index}, i, size(dis_p2g_all,2), t);
        else
            fprintf('the UDLF of %s has occur an error!\n',method{method_index});
            return;
        end

        reranking_result_matrix = load('output.txt');%%加载排序结果矩阵
        reranking_result = reranking_result_matrix(1,2:end);%%读取结果

        reranking_distance(:,i) = reranking_result';

    end

    save(['reranking_distance_',method{method_index}],'reranking_distance');
end


 %%画CMC曲线
for method_index = 1:length(method)
    
    load(['reranking_distance_',method{method_index}]);
    init_result = dis_p2g_all;
    ground_truth = repmat((1:size(init_result,2)),size(init_result,1),1);
    init_cmc = result_cmc(init_result,ground_truth);
    reranking_cmc = result_cmc(reranking_distance,ground_truth);
    figure(method_index);
    plot((1:316),init_cmc,'k');
    hold on
    h = plot((1:316),reranking_cmc,'b');
    title(method{method_index});
    saveas(h,[method{method_index},'.jpg']);
    
end

