function [score, icmr_time] = icmr17(dist, method, dataset_size)
% icmr17 7 method for reranking

method_list = {'CPRR','RLRECOM','RKGRAPH','CORGRAPH','RLSIM','RECKNNGRAPH','CONTEXTRR'};
assert(ismember(method, method_list));

%% config
test_dir = './method/baseline/icmr17/bin/';
exe_path = [test_dir 'udlf']; 

dist_dir = [test_dir 'dist.txt'];
output_dir = [test_dir 'output.txt'];
config_dir = [test_dir 'test_config.ini'];

%% initial

% write dist matrix into input.txt
% file = fopen(dist_dir, 'w');
% assert(file~=-1);
% 
% for j = 1:size(dist,1)
%     for k = 1:size(dist,2)
%         fprintf(file, '%f ', dist(j,k));
%     end
%     fprintf(file,'\n ');
% end
% fclose(file);

dlmwrite(dist_dir, dist, 'precision', '%5f', 'delimiter', '\t')
%% run
icmrtime  = tic;
[status, ~] = system([exe_path,' ',config_dir, ' ', method, ' ', num2str(dataset_size)]); 
icmr_time  = toc(icmrtime);

assert(status==0);

% collect re-ranking result (in dist format) into result
result = load(output_dir);
new_dist = result(1,2:end)';
score = exp(-new_dist);
