clear
load('pcm14_cuhk.mat');
% load('feature_cam1.mat');
% load('feature_cam2.mat');

gallery_num = 200;
query_num = 200;

feature = [mean(probe_feat,3)';mean(gallery_feat,3)'];
% feature = [feature_cam1;feature_cam2];

original_dist = MahDist(1, feature, feature); 
dist = re_ranking(original_dist,200,20,6,0.3);

final_dist = dist(1:query_num,query_num+1:end)';



init_result = MahDist(1,mean(probe_feat,3)',mean(gallery_feat,3)');
ground_truth = repmat((1:size(init_result,2)),size(init_result,1),1);
init_cmc = result_cmc(init_result',ground_truth);
init_auc = cmc2auc(init_cmc);

reranking_cmc(:,1) = result_cmc(final_dist',ground_truth);
reranking_auc(1) = cmc2auc(reranking_cmc(:,1));


dist2 = re_ranking(dist,200,20,6,0.3);

final_dist2 = dist2(1:query_num,query_num+1:end)';


reranking_cmc(:,2) = result_cmc(final_dist2',ground_truth);
reranking_auc(2) = cmc2auc(reranking_cmc(:,2));

save('pcm14_KRNN','reranking_cmc','reranking_auc','init_cmc','reranking_auc');
figure(1);
plot((1:200),init_cmc,'k');
hold on
plot((1:200),reranking_cmc(:,1),'b');
hold on
plot((1:200),reranking_cmc(:,2),'r');