function jaccard_dist = divide_and_fuse(probFea, galFea, k1, k2, iter, part_nums, lambda, burst)
    %% input featureï¼?num * ndim. 
    %% ndim denotes the dimension of the feature, 
    %% num denotes the size of the query/database set

    %% re-ranking parameter settings
    % k1 = 5;
    % k2 = 2;
    % iter = 2;
    % part_nums = [1, 2, 3, 4]; %1,8
    % lambda = 0.1; %fixed
    % burst = 0.5;  %fixed

prob_num = size(probFea, 2);
gal_num = size(galFea, 2);

%% iteration on different number of sub-features
for kk = 1:length(part_nums)

    part_num = part_nums(kk);
    part_length = ceil(size(galFea, 1) / part_num);
    
    %% define the vector for fusion
    V_fusion = zeros(gal_num+prob_num, gal_num, 'single');
    
    %% iteration on each sub-feature
    for jj = 1 : part_num

        %tic;
        %% select jj-th sub-feature
        start_idx = (jj-1) * part_length + 1;
        end_inx = min (jj * part_length, size(galFea, 1));

        part_galFea = galFea(start_idx : end_inx, :);
        part_probFea = probFea(start_idx : end_inx, :);
        
        %% compute original distance matrix
        original_dist = MahDist(1, [part_galFea part_probFea]', part_galFea');
        
        part_galFea = [];
        part_probFea = [];

        %% iterative SCA (except the most outer iteration)
        if iter >1
            for i = 1:iter-1
                temp_dist = SCA_reid_knnsimi(original_dist, k1, k2);
                original_dist = temp_dist*lambda + original_dist*(1-lambda);
            end
        end
        temp_dist = [];

        original_dist = original_dist./ repmat(max(original_dist, [], 2), 1, size(original_dist, 2));
        [~, initial_rank] = sort(original_dist, 2, 'ascend');

        test_num = size(original_dist,1);
        
        %% define the encoded vector of the jj-th sub-feature
        V = zeros(size(original_dist), 'single');

        [~, reverse_rank] = sort(initial_rank, 2, 'ascend');
        wseq = repmat(1:k1,k1,1);
        for ii = 1:test_num
            idx_now = initial_rank(ii, 1:k1);%ÅÅÇ°K1µÄË÷Òý
            idx_nb = reverse_rank(idx_now, idx_now);
            w = sum(1./(idx_nb.*wseq),2)';   %% k-NN similarity
            V(ii, idx_now) = w/sum(w);
        end
        reverse_rank = [];

        %% local query expansion
        if k2 ~=1
            V_qe = zeros(size(V), 'single' );
            for i = 1:size(V, 1)
                V_qe(i, :) = single(mean(V(initial_rank(i, 1:k2), :), 1));
            end
            V = V_qe;
            V_qe = [];
        end
        initial_rank = [];

        V_fusion = V_fusion + V.^burst ;

        %t = toc;
        %fprintf('i = %d / %d, (%d ~ %d)/%d, time %.2f s\n', jj, part_num, start_idx, end_inx, size(galFea, 1), t);

    end

    %% aggregation
    V_fusion = V_fusion / part_num;
    V_fusion = V_fusion.^(1/burst);
    %trick: L1-norm to accelerate the speed
    % V_fusion = bs_normalize_col(V_fusion',1)';

    %% Inverted Index
    invIndex = cell(gal_num, 1);
    for i = 1:gal_num
        invIndex{i} = find(V_fusion(1:gal_num , i) ~=0);
    end

    jaccard_dist = zeros(size(original_dist), 'single');

    for i = gal_num+1 : test_num 
        temp_min = zeros(1, gal_num, 'single');
        indNonZero = find( V_fusion( i, : ) ~= 0 );
        indImages = invIndex( indNonZero );
        for j = 1 : length( indNonZero )
            temp_min( 1, indImages{j} ) = temp_min( 1, indImages{j} )...
                + single( min( V_fusion(i, indNonZero(j)), V_fusion(indImages{j}, indNonZero(j)) ) )';
        end
        jaccard_dist(i, :) = bsxfun(@minus, 1, temp_min./(2 - temp_min)); 
    end

    jaccard_dist = jaccard_dist(gal_num+1 : test_num, : );
end 
end




function dist = SCA_reid_knnsimi(dist, k1, k2)
% Self-contex aggregation(SCA) is used to re-ranking the distance matrix to
% improve performance.
%It's the offline version. 
%Notes:
%   1.Ensure distpath load the varible "dist_pp"
%   2.Similar to TCA, just without two layer fusion. 
%Parameters:
%   1.k1 denotes the num of neighbors used in re-ranking, it has
%    great correlation with the average instance num of class
%   2.k2 denotes the num of neighbors used in neighbor augmentation
%     if k2 == 1, there will be no neighbor augmentation.here will be no 
%     neighbor augmentation. It may be 6~13
%Autors:Song Bai, zzchust in 2015/11/10

%(nd+nq * nd)
[q_num, d_num] = size(dist);
if q_num < d_num
    fprintf('n(row) should not less than n(col) !\n');
    dist = dist';
    temp = q_num;
    q_num = d_num;
    d_num = temp;
end
[~, idx] = sort(dist, 2, 'ascend');
f = zeros(size(dist), 'single');
dist = dist./ repmat(max(dist, [], 2), 1, size(dist, 2));% element/max

%% initiaze the contextual feature f
wseq = repmat(1:k1,k1,1);
[~, idx_rev] = sort(idx, 2, 'ascend');
for ii = 1:size(dist, 1)
    idx_now = idx(ii, 1:k1);
    idx_nb = idx_rev(idx_now, idx_now);
    w = sum(1./(idx_nb.*wseq),2)';
    f(ii, idx_now) = w/sum(w);
end

%% Neighbor augmentation
if k2 ~=1
    f1 = zeros(size(f), 'single');
    % fprintf('relevance feed back for k2 = %d...\n', k2);
    for ii = 1:size(f, 1)
        f1(ii, :) = single(mean(f(idx(ii, 1:k2), :)));
    end
    f = f1;
end
%% indexing
invIndex = cell(d_num, 1);
for ii = 1:d_num
    invIndex{ii} = find(f(1:d_num , ii) ~=0);
end
dist = zeros(q_num, d_num, 'single');
for i = 1:q_num
    temp_min = zeros(1, d_num, 'single');
    indNonZero = find( f( i, : ) ~= 0 );
    indImages = invIndex( indNonZero );
    for j = 1 : length( indNonZero )
        temp_min( 1, indImages{j} ) = temp_min( 1, indImages{j} )...
            + single( min( f(i, indNonZero(j)), f(indImages{j}, indNonZero(j)) ) )';
    end
    dist(i, :) = bsxfun(@minus, 1, temp_min./(2 - temp_min));
end

%If you find the NN is lower than the result of dist1 or dist2, the reason
%is through the re-ranking procedure the idx2(:,1) has been queued to tha
%last, you can force it to the front by a trick.
%% ensure NN
if 0
 for ii = 1:size(dist, 1)
     dist(ii, idx(ii, 2)) = 0.001;
 end
end
end