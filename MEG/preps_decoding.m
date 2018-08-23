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
if ~exist('subj',        'var'), subj         = 'pilot-005';                end
if ~exist('root_dir',    'var'), root_dir     = '/project/3011210.01/MEG/'; end
if ~exist('save_dir',    'var'), save_dir     = '/project/3011210.01/MEG/Classification'; end
if ~exist('classes',     'var'), classes      = {'ART', 'NN'};              end
if ~exist('classifier',  'var'), classifier   = 'preps_naivebayes';         end
if ~exist('timestep',    'var'), timestep     = 0.1;                        end
if ~exist('repeats',     'var'), repeats      = 100;                        end
if ~exist('folds',       'var'), folds        = 20;                         end
if ~exist('numfeat',     'var'), numfeat      = 250;                        end

%options to execute
if ~exist('dopca',       'var'), dopca        = true;                       end
if ~exist('donormal',    'var'), donormal     = true;                       end
if ~exist('doshuffle',   'var'), doshuffle    = false;                      end
if ~exist('compute_lc',  'var'), compute_lc   = false;                      end
if ~exist('compute_acc', 'var'), compute_acc  = false;                      end
if ~exist('dogeneral',   'var'), dogeneral    = false;                      end

if strcmp(subj,'pilot-002')
    warning('need to adjust trigger info')
end
%find which trigger number belongs to specified class;
for c = 1:length(classes)
    trig{c} = trigger{strcmp(pos,classes{c})};
end
trig = horzcat(trig{:});

if dogeneral
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

for c = 1:length(classes)
    indsel          = ismember(data.trialinfo(:,1),trigger{strcmp(pos,classes{c})});
    if dogeneral, indsel = indsel | ismember(data.trialinfo(:,1),testtrig(testlabel==c));end
    fprintf('selecting %d samples for class %d\n',sum(indsel),c)
    labels{c}       = ones(1,sum(indsel))*c;
    cfg             = [];
    cfg.trials      = indsel;
    datasel{c}      = ft_selectdata(cfg,data);
end
datasel = ft_appenddata([],datasel{:});
labels  = horzcat(labels{:})';

%check if equal amount of samples per class, if not resample
if ~exist('do_resample', 'var')
    if length(unique(sum(labels==labels'))) ~= 1
        do_resample = 1;
    else
        do_resample = 0;
    end
end

%% prepare data
%decompose matrix and only keep 60 components
if dopca
    cfg                   = [];
    cfg.demean            = 'yes';
    cfg.scale             =  0;
    cfg.method            = 'pca';
    cfg.numcomponent      = 60;
    datasel             = ft_componentanalysis(cfg, datasel);
end

cfg                   = [];
cfg.keeptrials        = 'yes';
cfg.vartrllength      = 1;
avg_data              = ft_timelockanalysis(cfg,datasel);

%% Configuration
cfg                   = [];
cfg.method            = 'crossvalidate';
cfg.mva               = classifier;
cfg.statistic         = {'accuracy'};
cfg.type              = 'nfold'; %'bloo' only with evenly distributed classes
cfg.nfolds            = folds;
cfg.resample          = do_resample;% default: false; resampling for 'loo' retuns zero
cfg.poolsigma         = 0;
cfg.numfeat           = numfeat;
cfg.design            = labels;

tsteps                = [avg_data.time(1)+timestep:timestep:avg_data.time(end)];

%% if compute accuracy
if compute_acc
    %% Loop over time & repeats
    %
    acc                   = zeros(length(tsteps),repeats);
    accshuf               = zeros(length(tsteps),repeats);
    
    for  t = 1:length(tsteps)
        cfgt            = [];
        cfgt.latency    = [tsteps(t)-timestep tsteps(t)];
        datatmp         = ft_selectdata(cfgt,avg_data);
        rng('default'); % ensure same 'random' behaviour for each time slice.
        for rep = 1:repeats
            if donormal
                out               = ft_timelockstatistics(cfg,datatmp);
                acc(t,rep)   = out.statistic.accuracy;
            end
            if doshuffle
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
                
                outshuf           = ft_timelockstatistics(cfgtmp,datatmp);
                accshuf(t,rep)    = outshuf.statistic.accuracy;
            end
        end
    end
    %%save results to file including timesteps
    cfg.timeinfo = tsteps;
    if donormal
        filename = fullfile(save_dir, subj, sprintf('classacc_%s_%dfolds_%dfeats_%s',subj,folds,numfeat,horzcat(classes{:})));
        save(filename, 'acc','cfg');
    end
    if doshuffle
        filename = fullfile(save_dir, subj, sprintf('classacc_%s_%dfolds_%dfeats_%s_shuf',subj,folds,numfeat,horzcat(classes{:})));
        save(filename, 'accshuf','cfg');
    end
end

%% if compute learning curve
if compute_lc
    if doshuffle
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
    filename = fullfile(save_dir, subj, sprintf('classlc_%s_%dfolds_%dfeats_%s',subj,folds,numfeat,horzcat(classes{:})));
    save(filename, 'acctest','acctrain','cfg');
end

%% If do generalize over time
if dogeneral
    
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
        
        if donormal
            out             = ft_timelockstatistics(cfgtmp,datatmp);
            acc(t)      = out.statistic.accuracy;
            acctrain(t) = out.trainacc.statistic.accuracy;
        end
        if doshuffle
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
    if donormal
        filename = fullfile(save_dir, subj, sprintf('classgeneral_%s_%dfeats_%sto%s',subj,numfeat,horzcat(classes{:}),horzcat(testpos{:})));
        save(filename, 'acc','acctrain','cfgtmp','prev');
    end
    if doshuffle
        filename = fullfile(save_dir, subj, sprintf('classgeneral_%s_%dfeats_%sto%s_shuf',subj,numfeat,horzcat(classes{:}),horzcat(testpos{:})));
        save(filename, 'accshuf','accshuftrain','cfgtmp','prev');
    end
end



