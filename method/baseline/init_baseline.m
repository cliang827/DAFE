baseline_para.method_set = {'mr', 'emr'};
baseline_para.tot_query_times = 3;

% parameters for manifold ranking (mr)
baseline_para.mr.k = 5;
baseline_para.mr.alpha = 0.99;

% parameters for efficient manifold ranking (emr)
baseline_para.emr.alpha = 0.99;
baseline_para.emr.p = 10;