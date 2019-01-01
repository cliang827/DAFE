function sim_score = dist2sim(dist_mat)

sigma = repmat(mean(dist_mat), size(dist_mat,1), 1);
feedback_score = exp(-1*dist_mat./sigma);
sim_score = normalization(feedback_score, [-1 1], 0, 'range-priority');


