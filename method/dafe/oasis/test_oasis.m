% test_robot_labeler

clear all
clc

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


    %% construct robot
    load_features_all; % load all features.
    t3 = tic;
    feedback_dist_set = cell(sys.setnum,1);
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

        parms.do_sym = 1;
        parms.use_matlab = 1;
        n = 1:4;
        model = oasis(feature(n,:), labels(n), parms);
        
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

%         dists = dists';
%         sigma = repmat(mean(dists), size(dists,1), 1);
%         feedback_score = exp(-1*dists./sigma);
%         feedback_score = normalization(feedback_score, [-1 1], 0, 'range-priority');
%         feedback_score_set{set,1} = feedback_score;
        feedback_dist_set{set,1} = dists';
    end
    time = toc(t3);
    fprintf(1, 'total time: %.0f sec.\n', time);

    save_dir = sprintf('./data/%s_gog_xqda.mat',databasename);
    save(save_dir, ...
        'testimagenames_set', 'testcamIDs_set', ...
        'probe_feat_set', 'gallery_feat_set', ...
        'g2p_dist_set', 'g2g_dist_set', ...
        'g2p_sim_set', 'g2g_sim_set', ...
        'feedback_dist_set', 'groundtruth_rank_set', '-v7.3');
end




 