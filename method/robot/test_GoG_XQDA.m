% test_robot_labeler

clear all
clc

% addpath('GOG/');
% addpath('GOG/mex');
% addpath('XQDA/');
% addpath('DataManage/');
% addpath('cliang_test/');
show_step = 50;

%% configuraiton of datasets. 
% 1 -- VIPeR,  2 -- CUHK01(M=1),  3 -- CUHK01(M=2),  4 -- PRID450s, 
% 5 -- GRID,  6 -- CUHK03(labeled),  7 -- CUHK03(detected)
for database_id = 1%[1 2 6 7 4 5]
    sys.database = database_id;
    
    if sys.database == 6 || sys.database == 7
        sys.setnum = 20;    % random division times
    else
        sys.setnum = 10;
    end
    set_database;

    % image size for resize
    H0 = 128; W0 = 48;

    %% configuration of features.  
    featuresetting = 1; % 1-- GOG_RGB, 2 -- GOG_Fusion
    parFea.featurenum = 4; 
    parFea.featureConf = cell( parFea.featurenum, 1);
    parFea.usefeature = zeros( parFea.featurenum, 1);

    switch featuresetting
        case 1 % GOG_RGB
            parFea.usefeature(1) = 1; % GOG_RGB
            parFea.usefeature(2) = 0; % GOG_Lab
            parFea.usefeature(3) = 0; % GOG_HSV
            parFea.usefeature(4) = 0; % GOG_nRnG
        case 2 % GOG_Fusion
            parFea.usefeature(1) = 1; % GOG_RGB
            parFea.usefeature(2) = 1; % GOG_Lab
            parFea.usefeature(3) = 1; % GOG_HSV
            parFea.usefeature(4) = 1; % GOG_nRnG
        otherwise
            fprintf('Undefined feature setting \n');
    end

    for f = 1:parFea.featurenum
        parFea.featureConf(f) = {set_default_parameter(f)};
    end


%     %% extract feature for all images.
%     fprintf('*** low level feature extraction *** \n');
%     fprintf('database = %s \n', databasename );
%     fprintf('+ extract feature for all images \n');
%     for f = 1:parFea.featurenum
%         if parFea.usefeature(f) == 1
%             param = parFea.featureConf{f};
%             fprintf('feature = %d [ %s ] \n', f, param.name);
%             feature_all = zeros( allimagenums, param.dimension );
% 
%             t0 = tic;
%             nchar = 0;
%             for imgind = 1:allimagenums
%                 if mod(imgind, show_step) == 0
%                     nchar = fprintf('imgind = %d / %d \n', imgind, allimagenums); 
%                 end
% 
%                 X = imread( strcat(datadirname, allimagenames{imgind})); % load image
%                 if size(X, 1) ~= H0 || size(X, 2) ~= W0; X = imresize(X, [H0 W0]); end % resize
% 
%                 feature_all(imgind, :) = GOG(X, param); % extract GOG
% 
%                 if mod(imgind, show_step) == show_step-1
%                     fprintf(1, repmat('\b', 1, nchar));
%                 end
%             end
%             fprintf(1, repmat('\b', 1, nchar));
%             fprintf(1, 'done!\n');
% 
%             feaTime = toc(t0);
%             meanTime = feaTime/allimagenums;
%             fprintf('mean feature extraction time %.3f seconds per image\n', meanTime);
% 
%             name1 = sprintf('feature_all_%s', param.name );
%             name = strcat( featuredirname, databasename, '_',  name1, '.mat');
%             fprintf('%s \n', name);
% 
%             save( name,  'feature_all', '-v7.3' );
%         end
%     end

    %% compute g2p_dist and g2g_dist
    t2 = tic;
    load_features_all; % load all features.
    tot = 2;
    
    g2p_dist_set = cell(sys.setnum,1);
    g2g_dist_set = cell(sys.setnum,1);
    probe_feat_set = cell(sys.setnum,1);
    gallery_feat_set = cell(sys.setnum,1);
    groundtruth_rank_set = cell(sys.setnum,1);

    for set=1:sys.setnum
        fprintf(1, 'set %02d : ', set);
        
        extract_feature_cell_from_all; % load test data
        conc_feature_cell; % feature concatenation
        
        % normalization
        temp = repmat(sqrt(diag(feature*feature')),1,size(feature,2));
        feature = feature./temp;

        camIDs = testcamIDs_set{set};
        probe_feat = feature(camIDs == 1, :);
        gallery_feat = feature(camIDs == 2, :);    
        probe_set_num = size(probe_feat,1);
        gallery_set_num = size(gallery_feat,1);

        % groundtruth_rank        
        labels = testlabels_set{set};
        probe_label = labels(camIDs == 1, :);
        gallery_label = labels(camIDs == 2, :);
        groundtruth_rank = zeros(gallery_set_num, probe_set_num);
        for i = 1:probe_set_num
            ix = find(gallery_label==probe_label(i));
            groundtruth_rank(:,i) = ix;
        end
        groundtruth_rank_set{set} = groundtruth_rank;

        % g2p_dist & g2g_dist
        fprintf(1, 'compute g2p_dist and g2g_dist...');
        nchar = 0;
        g2p_dist = zeros(gallery_set_num, probe_set_num);
        g2g_dist = zeros(gallery_set_num, gallery_set_num);
        for i = 1:gallery_set_num
            if mod(i, show_step) == 0
                nchar = fprintf('%.0f%% (%d/%d)\n', 100*i/gallery_set_num, i, gallery_set_num); 
            end

            diff_feat = repmat(gallery_feat(i,:), probe_set_num, 1)-probe_feat;
            g2p_dist(i,:) = sqrt(diag(diff_feat*diff_feat'))';

            diff_feat = repmat(gallery_feat(i,:), gallery_set_num, 1)-gallery_feat;
            g2g_dist(i,:) = sqrt(diag(diff_feat*diff_feat'))';

            if mod(i, show_step) == show_step-1
                fprintf(1, repmat('\b', 1, nchar));
            end
        end
        fprintf(1, repmat('\b', 1, nchar));
        fprintf(1, 'done!\n');

        probe_feat_set{set} = probe_feat;
        gallery_feat_set{set} = gallery_feat;

        g2p_dist_set{set} = g2p_dist;
        g2g_dist_set{set} = g2g_dist;  

    end
    time = toc(t2);
    fprintf(1, 'total time: %.0f sec.\n', time);

    %% construct robot
    t3 = tic;
    robot_dist_set = cell(sys.setnum,1);
    for set = 1:(sys.setnum)
        fprintf('----------------------------------------------------------------------------------------------------\n');
        fprintf('set = %d \n', set);
        fprintf('----------------------------------------------------------------------------------------------------\n');

        %% Training data
        tot = 1;
        extract_feature_cell_from_all;  % load training data
        apply_normalization; % feature normalization
        conc_feature_cell; % feature concatenation

        % train XQDA metric learning
        camIDs = traincamIDs_set{set};
        probX = feature(camIDs == 1, :);
        galX = feature(camIDs == 2, :);
        labels = trainlabels_set{set};
        probLabels = labels(camIDs == 1);
        galLabels = labels(camIDs == 2);

        [XQDAresult.W, XQDAresult.M, inCov, exCov] = XQDA(galX, probX, galLabels, probLabels);
        clear camIDs probX galX probX galLabels probLabels

        %% Test data
        tot = 2;
        extract_feature_cell_from_all; % load test data
        apply_normalization; % feature normalization
        conc_feature_cell; % feature concatenation

        % apply XQDA metric learning
        camIDs = testcamIDs_set{set};
        probX = feature(camIDs == 1, :);
        galX = feature(camIDs == 2, :);
        labels = testlabels_set{set};
        labelsPr = labels(camIDs == 1);
        labelsGa = labels(camIDs == 2);

        dists = MahDist(XQDAresult.M, galX * XQDAresult.W, probX * XQDAresult.W)';

        CMC = zeros( numel(labelsGa), 1);
        for p=1:numel(labelsPr)
            score = dists(p, :);
            [sortscore, ind] = sort(score, 'ascend');

            correctind = find( labelsGa(ind) == labelsPr(p));
            CMC(correctind:end) = CMC(correctind:end) + 1;
        end
        CMC = 100.*CMC/numel(labelsPr);
        CMCs(set, :) = CMC;

        fprintf(' Rank1,  Rank5, Rank10, Rank15, Rank20\n');
        fprintf('%5.2f%%, %5.2f%%, %5.2f%%, %5.2f%%, %5.2f%%\n\n', CMC([1,5,10,15,20]));
        clear camIDs probX galX probX galLabels probLabels options XQDAresult

        robot_dist_set{set,1} = dists';
    end
    time = toc(t3);
    fprintf(1, 'total time: %.0f sec.\n', time);

    save_dir = sprintf('./data/%s_gog_xqda.mat',databasename);
    save(save_dir, ...
        'testimagenames_set', 'testcamIDs_set', ...
        'probe_feat_set', 'gallery_feat_set', ...
        'robot_dist_set', 'groundtruth_rank_set', '-v7.3');
end




 