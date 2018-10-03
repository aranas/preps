%% Binary classifier on some MEG data file and corresponding labels
% Classifiers available are:
% preps_naivebayes
% {dml.standardizer dml.svm};
% {dml.standardizer dml.blogreg};

%% Set variables
pos = {'ART','NN','VVFIN','ADJA','APPR','NA','VA','Fill'};
trigger = {[111,114,121,124,211,214,221,224], %Determiner
    [112,115,122,125,212,215,222,225],        %Nouns
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
if ~exist('datasuffix',     'var'), datasuffix   = '';                         end
if ~exist('timestep',       'var'), timestep     = 0.1;                        end
cfgcv                   = [];%config for call to ft_timelockstatistics later on
cfgcv.method            = 'crossvalidate';
%options to execute
if ~exist('dow2v',          'var'), dow2v       = false;                       end
if ~exist('docateg',        'var'), docateg      = false;                      end
if ~exist('compute_acc',         'var'), compute_acc       = false;            end
if ~exist('compute_lc',          'var'), compute_lc        = false;            end
if ~exist('compute_general',     'var'), compute_general   = false;            end
if docateg
    %default parameters for doing binary classification
    if ~exist('doshuffle_strat','var'), doshuffle_strat = true;                end
    if ~exist('doshuffle_rand', 'var'), doshuffle_rand = false;                end
    if ~exist('classes',        'var'), classes      = {'VA', 'NA'};           end
    if ~exist('classifier',     'var'), classifier   = 'preps_naivebayes';     end
    if ~exist('folds',          'var'), folds        = 20;                     end
    if ~exist('numfeat',        'var'), numfeat      = 250;                    end
    cfgcv.numfeat           = numfeat;
    if ~exist('statfun',        'var'), statfun      = {'accuracy'};           end
    if compute_general
        %specify training and testing time window/trials when generalising
        %over time
        if ~exist('trainwindow','var'), trainwindow  = [0.3 0.4];              end
        if ~exist('testtrig',   'var'), testtrig     = horzcat(trigger{6:7});  end
        fprintf(strcat('generalising to events ',repmat(' %d',size(testtrig)),'\n'),testtrig)
    end
elseif dow2v
    %default parameters for doing regression on continuous variable (ie
    %word embedding)
    if ~exist('doshuffle_strat','var'), doshuffle_strat = false;               end
    if ~exist('doshuffle_rand', 'var'), doshuffle_rand = true;                 end
    if ~exist('classes',        'var'), classes      = {'NN','VVFIN','ADJA','NA','VA'};           end
    if ~exist('classifier',     'var'), classifier   = 'ridgeregression_sa';   end
    if ~exist('statfun',        'var'), statfun      = {'eval_correlation'};   end
    if ~exist('lambda',         'var'), lambda      = [];                       end
    if ~exist('lambdaeval',     'var'), lambdaeval      = 'mse';              end
    cfgcv.lambdaeval = lambdaeval;
    cfgcv.lambda = lambda;
    %folds will be set according to number of trials
    suffix = '_w2v';
    do_resample = 0;
end
if ~compute_lc
    %compute 100 models for different folding schemes
    if ~exist('repeats',        'var'), repeats      = 100;                end
end
if ~exist('compute_tuda',        'var')
    compute_tuda      = false;
    clear numfeat
end

%options to prepare data
if ~exist('dopca',          'var'), dopca        = true;                       end
if ~exist('clean_muscle',   'var'), clean_muscle = false;                      end
if ~exist('dopretest100',   'var'), dopretest100 = false;                      end
if ~exist('doposttest',     'var'), doposttest   = false;                      end


if strcmp(subj,'pilot-002')
    warning('need to adjust trigger info')
end
%% Load & select data
load(strcat(root_dir,subj,'_dataclean',datasuffix,'.mat'))
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

if clean_muscle
    load(fullfile(root_dir,strcat(subj,'_muscle')));
    cfg         = [];
    cfg.trials  = ~ismember(data.trialinfo(:,3),noisy_trials(:,end));
    data        = ft_selectdata(cfg,data);
end

%find which trigger number belongs to specified class;
trig = cell(length(classes),1);
for c = 1:length(classes)
    trig{c} = trigger{strcmp(pos,classes{c})};
end
trig = horzcat(trig{:});
if compute_general
    if ~any(ismember(testtrig,trig)) && ~exist('testlabel','var')
        error('labels for test samples need to be specified')
    else
        fprintf(strcat('labeling events as ',repmat(' %d',size(testlabel)),'\n'),testlabel)
    end
end
% select all trials belonging to specified class and construct
% class-specific labels
labels = cell(length(classes),1);
datasel = cell(length(classes),1);
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
    labels = zeros(size(datasel.trial,2),300);
    indsel = false(1,size(datasel.trial,2));
    for smp = 1:size(datasel.trial,2)
        id          = datasel.trialinfo(smp,2);
        trig        = num2str(datasel.trialinfo(smp,1));
        w2vtmp      = stimuli(id).words(str2num(trig(end))).w2v;
        if size(w2vtmp,2)~=0
            labels(smp,:) = stimuli(id).words(str2num(trig(end))).w2v;
            indsel(smp) = 1;
        end
    end
    labels            = labels(indsel,:);
    cfg             = [];
    cfg.trials      = indsel;
    datasel         = ft_selectdata(cfg,datasel);
    
    % set some variables
    div                = divisors(size(labels,1));
    div                = div(mod(div,2)==0);
    folds              = size(labels,1)/div(6);
end


% check if equal amount of samples per class, if not resample
if ~exist('do_resample', 'var')
    if length(unique(sum(labels==labels'))) ~= 1
        do_resample = 1;
    else
        do_resample = 0;
    end
end

%clear some memory
clear data
%% Continue setting parameters in Configuration
cfgcv.mva               = classifier;
cfgcv.statistic         = statfun;
cfgcv.type              = 'nfold'; %'bloo' only with evenly distributed classes
cfgcv.nfolds            = folds;
cfgcv.poolsigma         = 0;
cfgcv.design            = labels;
cfgcv.resample          = do_resample;% default: false; resampling for 'loo' retuns zero
if size(timestep,2)>1
    tsteps = timestep;
    timestep     = 0.1;   
else
tsteps                = [datasel.time{1}(1)+timestep:timestep:datasel.time{1}(end)];
end
%% prepare data
%decompose matrix and only keep 60 components
if dopca
    cfgtmp                   = [];
    cfgtmp.demean            = 'yes';
    cfgtmp.scale             =  0;
    cfgtmp.method            = 'pca';
    cfgtmp.numcomponent      = 60;
    datasel                  = ft_componentanalysis(cfgtmp, datasel);
end

cfgtmp                   = [];
cfgtmp.keeptrials        = 'yes';
cfgtmp.vartrllength      = 1;
datasel                  = ft_timelockanalysis(cfgtmp,datasel);

%% if compute accuracy
if compute_acc
    %% Loop over time & repeats
    %
    
    for  t = 1:length(tsteps)
        cfgt            = [];
        cfgt.latency    = [tsteps(t)-timestep tsteps(t)];
        datatmp         = ft_selectdata(cfgt,datasel);
        rng('default'); % ensure same 'random' behaviour for each time slice.
        for rep = 1:repeats
            out               = ft_timelockstatistics(cfgcv,datatmp);
            if docateg
                acc(t,rep)   = out.statistic.accuracy;
            elseif dow2v
                acc.acc(t,rep)   = out.statistic{1};
                acc.model(t,rep,:,:) = out.model;
                acc.out{t,rep} = out.out;
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
            end
            if doshuffle_rand
                labels_perm       = labels(randperm(size(labels,1)),:);
            end
            if doshuffle_strat || doshuffle_rand
                cfgtmp            = cfgcv;
                cfgtmp.design     = labels_perm;
                outshuf           = ft_timelockstatistics(cfgtmp,datatmp);
                if docateg
                    accshuf(t,rep)   = outshuf.statistic.accuracy;
                elseif dow2v
                    cfg.lambda   = lambda;
                    accshuf(t,rep)   = outshuf.statistic{1};
                    model(t,rep,:,:) = out.model;
                end
            end
        end
    end
    %%save results to file including timesteps
    cfg.timeinfo = tsteps;
    if docateg
        filename = fullfile(save_dir, subj, sprintf('classacc_%s%s_%dfolds_%dfeats_%s%s',subj,datasuffix,folds,numfeat,horzcat(classes{:}),suffix));
        save(filename, 'acc','cfgcv');
    end
    if doshuffle_strat
        filename = fullfile(save_dir, subj, sprintf('classacc_%s%s_%dfolds_%dfeats_%s%s_shuf',subj,datasuffix,folds,numfeat,horzcat(classes{:}),suffix));
        save(filename, 'accshuf','cfgtmp');
    end
    if dow2v
        filename = fullfile(save_dir, subj, sprintf('classacc_%s%s_%dfolds_lambda%d_%s%s',subj,datasuffix,folds,cfgcv.lambda,horzcat(classes{:}),suffix));
        save(filename, 'acc','cfgcv');
    end
    if doshuffle_rand
        filename = fullfile(save_dir, subj, sprintf('classacc_%s%s_%dfolds_lambda%d_%s%s_shuf',subj,datasuffix,folds,cfgcv.lambda,horzcat(classes{:}),suffix));
        save(filename, 'accshuf','cfgtmp');
    end
end

%% if compute learning curve
if compute_lc
    N         = length(labels);
    grid_smp  = [2:4:N-N/20];
    
    %% Loop over time & repeats
    acctest                   = zeros(length(tsteps),length(grid_smp));
    acctrain                   = zeros(length(tsteps),length(grid_smp));
    for  t = 1:length(tsteps)
        cfgt            = [];
        cfgt.latency    = [tsteps(t)-timestep tsteps(t)];
        datatmp         = ft_selectdata(cfgt,datasel);
        rng('default'); % ensure same 'random' behaviour for each time slice.
        parfor nsmp = 1:size(grid_smp,2)
            cfgtmp             = cfgcv; %need to assign variable within parfor loop
            cfgtmp.max_smp     = grid_smp(nsmp);
            out                = ft_timelockstatistics(cfgtmp,datatmp);
            acctest(t,nsmp)    = out.statistic.accuracy;
            acctrain(t,nsmp)   = out.trainacc.statistic.accuracy;
        end
    end
    %%save results to file including timesteps
    cfg.timeinfo = tsteps;
    filename = fullfile(save_dir, subj, sprintf('classlc_%s%s_%dfolds_%dfeats_%s%s',subj,datasuffix,folds,numfeat,horzcat(classes{:}),suffix));
    save(filename, 'acctest','acctrain','cfgtmp');
end

%% if do generalize over time
if compute_general
    
    %remember order of trials for appenddata
    datasel.trialinfo(:,end+1) = [1:length(datasel.trialinfo)]';
    %select time window to be trained on
    cfgt             = [];
    cfgt.trials       = ~ismember(datasel.trialinfo(:,1),testtrig);
    cfgt.latency     = trainwindow;
    data_train       = ft_selectdata(cfgt,datasel);
    
    for t = 1:length(tsteps) % for each time window to be tested on
        cfgt              = [];
        cfgt.trials       = ismember(datasel.trialinfo(:,1),testtrig);
        cfgt.latency      = [tsteps(t)-timestep tsteps(t)];
        data_test        = ft_selectdata(cfgt,datasel);
        data_test.time   = data_train.time;
        
        cfgt              = [];
        datatmp           = ft_appenddata(cfgt,data_train,data_test);
        %recover original trial order
        [~,idx]           = sort(datatmp.trialinfo(:,end));
        datatmp.trial     = datatmp.trial(idx);
        datatmp.trialinfo = datatmp.trialinfo(idx,:);
        %set parameters
        cfgtmp            = cfgcv;
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
        filename = fullfile(save_dir, subj, sprintf('classgeneral_%s%s_%dfeats_%sto%s%s',subj,datasuffix,numfeat,horzcat(classes{:}),horzcat(testpos{:}),suffix));
        save(filename, 'acc','acctrain','cfgtmp','prev');
    end
    if doshuffle_strat
        filename = fullfile(save_dir, subj, sprintf('classgeneral_%s%s_%dfeats_%sto%s%s_shuf',subj,datasuffix,numfeat,horzcat(classes{:}),horzcat(testpos{:}),suffix));
        save(filename, 'accshuf','accshuftrain','cfgtmp','prev');
    end
end

%% if do Hidden markov model as described in Vidaurre et al. 2018
if compute_tuda
    cfg.constant = 0;
    K = 4;
    
   
    [N p ttrial] = size(datasel.trial);
    datatmp = reshape(permute(datasel.trial,[3 1 2]),[ttrial*N p]);
    datatmp(any(isnan(datatmp),2),:) = [];
    T =repmat(length(datasel.time),[N 1]);
    
    options.K = K;
    options.parallel_trials = 1;
    [tuda,Gamma] = tudatrain (datatmp,labels,T,options);
    
    cfgcv.Gamma = Gamma; 
    out               = ft_timelockstatistics(cfgcv,datasel);
    
end



