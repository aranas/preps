%% Binary classifier on some MEG data file and corresponding labels
% Classifiers available are:
% preps_naivebayes
% {dml.standardizer dml.svm};
% {dml.standardizer dml.blogreg};

%% Set variables
pos = {'ART','NN','VVFIN','ADJA','APPR','NA','VA','Fill'};
trigger = {[111,114,121,124,211,214,221,224], %Determiner
    [112,115,122,125,212,215,222,225,119,129,219,229],        %Nouns
    [113,123,213,223],                        %Verbs
    [118,128,218,228],                        %Adjectives
    [116,126,216,226],                        %Preposition
    [219,229],                                %last word Noun attached
    [119,129],                                %last word Verb attached
    [30:39]};                                 %all words in filler sentences
%parameters
if ~exist('subj',           'var'), subj         = 'pilot-005';                end
if ~exist('root_dir',       'var'), root_dir     = '/project/3011210.01/MEG/'; end
if ~exist('save_dir',       'var'), save_dir     = '/project/3011210.01/MEG/Classification'; end
if ~exist('suffix',         'var'), suffix       = '';                         end %might change later in code depending on selected options
if ~exist('classes',        'var'), classes      = {'VA', 'NA'};              end
if ~exist('classifier',     'var'), classifier   = 'preps_naivebayes';         end
if ~exist('timestep',       'var'), timestep     = 0.1;                        end
if ~exist('repeats',        'var'), repeats      = 100;                        end
if ~exist('folds',          'var'), folds        = 20;                         end
if ~exist('numfeat',        'var'), numfeat      = 250;                        end
if ~exist('statfun',        'var'), statfun      = {'accuracy'};               end

%options to prepare data
if ~exist('dopca',          'var'), dopca        = true;                       end
if ~exist('docateg',        'var'), docateg      = false;                      end
if ~exist('dow2v',          'var'), dow2v       = false;                       end
if ~exist('doshuffle_rand', 'var'), doshuffle_rand = false;                    end
if ~exist('doshuffle_strat','var'), doshuffle_strat = false;                   end
if ~exist('dopretest100',   'var'), dopretest100 = false;                      end
if ~exist('doposttest',     'var'), doposttest   = false;                      end
%options to execute
if ~exist('compute_lc',          'var'), compute_lc        = false;                      end
if ~exist('compute_acc',         'var'), compute_acc       = false;                      end
if ~exist('compute_general',     'var'), compute_general   = false;                      end
if ~exist('compute_tuda',        'var'), compute_tuda      = false;                      end


if strcmp(subj,'pilot-002')
    warning('need to adjust trigger info')
end
%find which trigger number belongs to specified class;
for c = 1:length(classes)
    trig{c} = trigger{strcmp(pos,classes{c})};
end
trig = horzcat(trig{:});

if compute_general
    if ~exist('trainwindow','var'), trainwindow  = [0.3 0.4];               end
    if ~exist('testtrig',   'var'), testtrig     = horzcat(trigger{6:7});   end
    fprintf(strcat('generalising to events ',repmat(' %d',size(testtrig)),'\n'),testtrig)
    
    if ~any(ismember(testtrig,trig)) && ~exist('testlabel','var')
        error('labels for test samples need to be specified')
    else
        fprintf(strcat('labeling events as ',repmat(' %d',size(testlabel)),'\n'),testlabel)
    end
end
%% Load & select data
load(strcat(root_dir,subj,'_dataclean.mat'))
clear badcomp compds

% if dopretest100 only keep trials that scored high on pre-test
if dopretest100
    load preps_stimuli
    accs        = [stimuli(1:200).acc];
    idsel       = find(accs>=0.9);
    cfg         = [];
    cfg.trials  = ismember(data.trialinfo(:,2),idsel);
    data        = ft_selectdata(cfg,data);
    suffix      = '_pretest100';
end

% select all trials belonging to specified class and construct
% class-specific labels
for c = 1:length(classes)
    indsel          = ismember(data.trialinfo(:,1),trigger{strcmp(pos,classes{c})});
    if compute_general, indsel = indsel | ismember(data.trialinfo(:,1),testtrig(testlabel==c));end
    fprintf('selecting %d samples for class %d\n',sum(indsel),c)
    labels{c}       = ones(1,sum(indsel))*c;
    cfg             = [];
    cfg.trials      = indsel;
    datasel{c}      = ft_selectdata(cfg,data);
end
datasel = ft_appenddata([],datasel{:});
labels  = horzcat(labels{:})';

% check if equal amount of samples per class, if not resample
if ~exist('do_resample', 'var')
    if length(unique(sum(labels==labels'))) ~= 1
        do_resample = 1;
    else
        do_resample = 0;
    end
end

% if doposttest relabel samples according to individual post tests;
if doposttest
    %relabel samples according to post-test
    load preps_stimuli
    subjn = str2num(subj(end-2:end));
    for smp = 1:size(datasel.trial,2)
        id = datasel.trialinfo(smp,2);
        if strcmp(stimuli(id).post_test(subjn).attachment,'Nomen')
            labels(smp) = 1;
        else
            labels(smp) = 2;
        end
    end
    suffix      = '_posttest';
end

if dow2v
    load preps_stimuli
    
    % replace categorical labels with w2v info
    feat = zeros(size(datasel.trial,2),300);
    indsel = false(1,size(datasel.trial,2));
    for smp = 1:size(datasel.trial,2)
        id          = datasel.trialinfo(smp,2);
        trig        = num2str(datasel.trialinfo(smp,1));
        w2vtmp      = stimuli(id).words(str2num(trig(end))).w2v;
        if size(w2vtmp,2)~=0
            feat(smp,:) = stimuli(id).words(str2num(trig(end))).w2v;
            indsel(smp) = 1;
        end
    end
    feat            = feat(indsel,:);
    cfg             = [];
    cfg.trials      = indsel;
    datasel         = ft_selectdata(cfg,datasel);
    
    % overwrite some variables
    labels             = feat;
    classifier         = 'ridgeregression_sa';
    statfun            = {'eval_correlation'};
    div                = divisors(size(feat,1));
    folds              = size(feat,1)/div(3);
    do_resample        = 0;
    %cfg.lambda         = 7.9207e+06;
    suffix = '_w2v';
end

%% Configuration
cfg                   = [];
cfg.method            = 'crossvalidate';
cfg.mva               = classifier;
cfg.statistic         = statfun;
cfg.type              = 'nfold'; %'bloo' only with evenly distributed classes
cfg.nfolds            = folds;
cfg.resample          = do_resample;% default: false; resampling for 'loo' retuns zero
cfg.poolsigma         = 0;
cfg.numfeat           = numfeat;
cfg.design            = labels;

tsteps                = [datasel.time{1}(1)+timestep:timestep:datasel.time{1}(end)];
if strcmp(tmpstep,'all')
    tsteps            = datasel.time{1};
end
%% prepare data
%decompose matrix and only keep 60 components
if dopca
    cfgtmp                   = [];
    cfgtmp.demean            = 'yes';
    cfgtmp.scale             =  0;
    cfgtmp.method            = 'pca';
    cfgtmp.numcomponent      = 60;
    datasel             = ft_componentanalysis(cfgtmp, datasel);
end

cfgtmp                   = [];
cfgtmp.keeptrials        = 'yes';
cfgtmp.vartrllength      = 1;
avg_data              = ft_timelockanalysis(cfgtmp,datasel);

%% if compute accuracy
if compute_acc
    %% Loop over time & repeats
    %
    acc                   = zeros(length(tsteps),repeats);
    accshuf               = zeros(length(tsteps),repeats);
    
    for  t = 1:length(tsteps)
        cfgt            = [];
        cfgt.latency    = [tsteps(t)-timestep tsteps(t)];
        if strcmp(timestep,'all'), cfgt.latency = tsteps(t);  end
        datatmp         = ft_selectdata(cfgt,avg_data);
        rng('default'); % ensure same 'random' behaviour for each time slice.
        for rep = 1:repeats
            out               = ft_timelockstatistics(cfg,datatmp);
            if docateg
                acc(t,rep)   = out.statistic.accuracy;
            elseif dow2v
                acc(t,rep)   = out.statistic{1};
                model(t,rep,:,:) = out.model;
            end
            if doshuffle_strat
                %permute labels
                indx_1            = find(labels==1);
                indx_2            = find(labels==2);
                indx_1            = indx_1(randperm(size(indx_1,1)));
                indx_2            = indx_2(randperm(size(indx_2,1)));
                
                labels_perm       = labels;
                labels_perm(indx_1(1:(length(indx_1)/2))) = 2;
                labels_perm(indx_2(1:(length(indx_2)/2))) = 1;
                cfgtmp            = cfg;
                cfgtmp.design     = labels_perm;
            end
            if doshuffle_rand
                labels_perm       = labels(randperm(size(labels,1)),:);
            end
            if doshuffle_strat || doshuffle_rand
                cfgtmp            = cfg;
                cfgtmp.design     = labels_perm;
                outshuf           = ft_timelockstatistics(cfgtmp,datatmp);
                if docateg
                    acc(t,rep)   = outshuf.statistic.accuracy;
                elseif dow2v
                    acc(t,rep)   = outshuf.statistic{1};
                    model(t,rep,:,:) = out.model;
                end
            end
        end
    end
    %%save results to file including timesteps
    cfg.timeinfo = tsteps;
    if docateg
        filename = fullfile(save_dir, subj, sprintf('classacc_%s_%dfolds_%dfeats_%s%s',subj,folds,numfeat,horzcat(classes{:}),suffix));
        save(filename, 'acc','cfg');
    end
    if doshuffle_strat || doshuffle_rand
        filename = fullfile(save_dir, subj, sprintf('classacc_%s_%dfolds_%dfeats_%s%s_shuf',subj,folds,numfeat,horzcat(classes{:}),suffix));
        save(filename, 'accshuf','cfg');
    end
end

%% if compute learning curve
if compute_lc
    if doshuffle_strat
        warning ('cannot compute both shuffles/repeats & learning curve')
    end
    repeats   = 1;
    m         = length(labels);
    grid_smp  = [2:4:m-m/20];
    
    %% Loop over time & repeats
    acc                   = zeros(length(tsteps),length(grid_smp));
    
    for  t = 1:length(tsteps)
        cfgt            = [];
        cfgt.latency    = [tsteps(t)-timestep tsteps(t)];
        datatmp         = ft_selectdata(cfgt,avg_data);
        rng('default'); % ensure same 'random' behaviour for each time slice.
        parfor nsmp = 1:size(grid_smp,2)
            cfgtmp             = cfg; %need to assign variable within parfor loop
            cfgtmp.max_smp     = grid_smp(nsmp);
            out                = ft_timelockstatistics(cfgtmp,datatmp);
            acctest(t,nsmp)    = out.statistic.accuracy;
            acctrain(t,nsmp)   = out.trainacc.statistic.accuracy;
        end
    end
    %%save results to file including timesteps
    cfg.timeinfo = tsteps;
    filename = fullfile(save_dir, subj, sprintf('classlc_%s_%dfolds_%dfeats_%s%s',subj,folds,numfeat,horzcat(classes{:}),suffix));
    save(filename, 'acctest','acctrain','cfg');
end

%% if do generalize over time
if compute_general
    
    %remember order of trials for appenddata
    avg_data.trialinfo(:,end+1) = [1:length(avg_data.trialinfo)]';
    %select time window to be trained on
    cfgt             = [];
    cfgt.trials       = ~ismember(avg_data.trialinfo(:,1),testtrig);
    cfgt.latency     = trainwindow;
    data_train       = ft_selectdata(cfgt,avg_data);
    
    for t = 1:length(tsteps) % for each time window to be tested on
        cfgt              = [];
        cfgt.trials       = ismember(avg_data.trialinfo(:,1),testtrig);
        cfgt.latency      = [tsteps(t)-timestep tsteps(t)];
        data_test        = ft_selectdata(cfgt,avg_data);
        data_test.time   = data_train.time;
        
        cfgt              = [];
        datatmp           = ft_appenddata(cfgt,data_train,data_test);
        %recover original trial order
        [~,idx]           = sort(datatmp.trialinfo(:,end));
        datatmp.trial     = datatmp.trial(idx);
        datatmp.trialinfo = datatmp.trialinfo(idx,:);
        %set parameters
        cfgtmp            = cfg;
        cfgtmp.type       = 'split';
        cfgtmp.max_smp    = sum(~ismember(datatmp.trialinfo(:,1),testtrig));
        cfgtmp.testfolds  = {find(ismember(datatmp.trialinfo(:,1),testtrig))};
        
        if docateg
            out             = ft_timelockstatistics(cfgtmp,datatmp);
            acc(t)      = out.statistic.accuracy;
            acctrain(t) = out.trainacc.statistic.accuracy;
        end
        if doshuffle_strat
            parfor rep = 1:repeats
                indx_1            = find(labels==1);
                indx_2            = find(labels==2);
                indx_1            = indx_1(randperm(size(indx_1,1)));
                indx_2            = indx_2(randperm(size(indx_2,1)));
                
                labels_perm       = labels;
                labels_perm(indx_1(1:(length(indx_1)/2))) = 2;
                labels_perm(indx_2(1:(length(indx_2)/2))) = 1;
                cfginner            = cfgtmp;
                cfginner.design     = labels_perm;
                
                outshuf             = ft_timelockstatistics(cfginner,datatmp);
                accshuf(t,rep)      = outshuf.statistic.accuracy;
                accshuftrain(t,rep) = outshuf.trainacc.statistic.accuracy;
            end
        end
    end
    %%save results to file including timesteps
    cfgtmp.timeinfo    = tsteps;
    cfgtmp.trainwindow = trainwindow;
    cfgtmp.testtrigger = testtrig;
    cfgtmp.testlabel   = testlabel;
    prev = datatmp.cfg.previous;
    testpos = pos(cellfun(@(x) any(ismember(x,testtrig)),trigger));
    if docateg
        filename = fullfile(save_dir, subj, sprintf('classgeneral_%s_%dfeats_%sto%s%s',subj,numfeat,horzcat(classes{:}),horzcat(testpos{:}),suffix));
        save(filename, 'acc','acctrain','cfgtmp','prev');
    end
    if doshuffle_strat
        filename = fullfile(save_dir, subj, sprintf('classgeneral_%s_%dfeats_%sto%s%s_shuf',subj,numfeat,horzcat(classes{:}),horzcat(testpos{:}),suffix));
        save(filename, 'accshuf','accshuftrain','cfgtmp','prev');
    end
end

%% if do Hidden markov model as described in Vidaurre et al. 2018
if compute_tuda
    modelfile = fullfile(save_dir,subj,strcat('classacc_',subj,'_20folds_250feats_NNVVFINADJA_w2v.mat'));
    if ~exist(modelfile, 'file')
        qsubfeval('preps_execute_pipeline','preps_decoding',{'subj',subj},{'classes',classes},{'compute_acc',1},{'repeats',1},{'dow2v',1},{'timestep','all'},'memreq',10*1024^3,'timreq',3*60*60,'batchid',sprintf('preps_decoding_w2v_allt%s',subj));
        error('models have not been computed yet, job has been deployed')
    end
    
    load(modelfile)

    %compute error of each model e(t,j) = sum(sqrt(X(j)v(t) - Y(j)));
    
    %compute divergence between model v(t) and v(j);
    
    %weights = group T decoding models into K clusters (hierarchical clustering)
    
    %refine weights by expectations maximisation
    
    %get rid of synchrony across trials using bayesian appraoch (hmm-mar
    %toolbox)
end



