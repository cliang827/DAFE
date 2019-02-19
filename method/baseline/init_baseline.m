load(ctrl_para.dir_info.data_file_dir, 'probe_feat_set', 'gallery_feat_set');

baseline_para.method_set = {'mr', 'emr'};
baseline_para.method_set_icmr17 = {'CPRR','RLRECOM','RKGRAPH',...
    'CORGRAPH','RLSIM','RECKNNGRAPH','CONTEXTRR'};
baseline_para.method_set_other = {'krnn'};

baseline_para.tot_query_times = ctrl_para.exp.tot_query_times;

% parameters for manifold ranking (mr)
baseline_para.mr.k = 10;
baseline_para.mr.alpha = 0.99;

% parameters for efficient manifold ranking (emr)
baseline_para.emr.alpha = 0.99;
baseline_para.emr.p = 10;

%parameters for k-reciprocal reranking (krnn)
baseline_para.krnn.k1 = 20;
baseline_para.krnn.k2 = 6;
baseline_para.krnn.lambda = 0.3;

baseline_para_set = cell(1,para_test_num);
for i=1:para_test_num
    t = para_test_set{i,7};
    baseline_para.probe_feat = probe_feat_set{t};
    baseline_para.gallery_feat = gallery_feat_set{t};
    
    if debug_flag
        baseline_para.probe_feat = baseline_para.probe_feat(1:curr_dataset.probe_set_num,:);
    end
    
	baseline_para_set{i} = baseline_para;
end