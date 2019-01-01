function M = similarity_learning(M, pairwise_ix_set, feature_set, method)

n = size(pairwise_ix_set,1);
p = feature_set(end,:);
C = 1;

switch method
    case 'oasis'
        for set=1:n
            pos_ix = pairwise_ix_set(set,1);
            neg_ix = pairwise_ix_set(set,2);
            score_diff = pairwise_ix_set(set,3);

            samples_delta = [+1, -1] * feature_set([pos_ix, neg_ix],:);
            loss = score_diff - p*M*samples_delta';
            if loss>0
                grad_W = p'*samples_delta;
                norm_grad_W = sum(p .* p) * sum(samples_delta .*samples_delta);


                % constraint on the maximal update step size
                tau_val = loss/norm_grad_W; %loss / (V*V');
                tau = min(C, tau_val); 

                M = M + tau*grad_W;
            end
            loss = score_diff - p*M*samples_delta';
        end
end

M = 0.5*(M+M');