%% Binary classifier on some MEG data file and corresponding labels
% Classifiers available are:
% preps_naivebayes
% {dml.standardizer dml.svm};
% {dml.standardizer dml.blogreg};

%% Default parameters
if ~exist('subj',        'var'), subj         = 'pilot-005';                end
if ~exist('root_dir',    'var'), root_dir     = '/project/3011210.01/MEG/'; end
if ~exist('classes',     'var'), classes      = {'ART', 'NN'};              end
if ~exist('classifier',  'var'), classifier   = 'preps_naivebayes';         end
if ~exist('timestep',    'var'), timestep     = 0.1;                        end
if ~exist('repeats',     'var'), repeats      = 50;                         end
if ~exist('folds',       'var'), folds        = 20;                         end
if ~exist('numfeat',     'var'), folds        = 250;                        end
if ~exist('dopca',       'var'), dopca        = true;                       end
if ~exist('doshuffle',   'var'), doshuffle    = true;                       end
if ~exist('doqsub',      'var'), doqsub       = true;                       end

pos = {'ART','NN','VVFIN','ADJA','APPR','NA','VA','Fill'};
trigger = {[111,114,121,124,211,214,221,224], %Determiner
    [112,115,122,125,212,215,222,225],        %Nouns
    [113,123,213,223],                        %Verbs
    [118,128,218,228],                        %Adjectives
    [116,126,216,226],                        %Preposition
    [219,229],                                %last word Noun attached
    [119,129],                                %last word Verb attached
    [30:39]};                                 %all words in filler sentences  
if strcmp(subj,'pilot-002')
    warning('need to adjust trigger info')
end
%% Load & select data
load(strcat(root_dir,subj,'_dataclean.mat'))
clear badcomp compds

for c = 1:length(classes)
    indsel          = ismember(data.trialinfo(:,1),trigger{strcmp(pos,classes{c})});
    fprintf('selecting %d samples for class %d\n',sum(indsel),c)
    labels{c}       = ones(1,sum(indsel))*c;
    cfg             = [];
    cfg.trials      = indsel;
    datasel{c}      = ft_selectdata(cfg,data);
end
datasel = ft_appenddata([],datasel{:});
labels  = horzcat(labels{:});
%% prepare data
%decompose matrix and only keep 60 components
if dopca
    cfg                   = [];
    cfg.demean            = 'yes';
    cfg.scale             =  0;
    cfg.method            = 'pca';
    cfg.numcomponent      = 60;
    data_comp             = ft_componentanalysis(cfg, datasel);
end

cfg                   = [];
cfg.keeptrials        = 'yes';
avg_data              = ft_timelockanalysis(cfg,data_comp);

%% Configuration
cfg                   = [];
cfg.method            = 'crossvalidate';
cfg.mva               = classifier;
cfg.statistic         = {'accuracy'};
cfg.type              = 'nfold'; %'bloo' only with evenly distributed classes
cfg.nfolds            = folds;
cfg.resample          = 1;% default: false; resampling for 'loo' retuns zero
cfg.design            = labels;
cfg.poolsigma         = 0;
%cfg.max_smp           = 0;
cfg.numfeat           = 250; 
%cfg.randomseed        = 5;

%% Loop over time & repeat
%
tsteps                = [timestep:timestep:avg_data.time(end)];
jobid                 = cell(length(tsteps),repeat);
acc                   = zeros(length(tsteps),repeat);

for  t = 1:length(tsteps)
      cfgt            = [];
      cfgt.latency    = [tsteps(t)-timestep tsteps(t)];
      datatmp         = ft_selectdata(cfgt,avg_data);
      
  for rep = 1:repeats
    out             = ft_timelockstatistics(cfg,datatmp);
    %permute labels
    indx_1            = find(labels==1);
    indx_2            = find(labels==2);
    indx_1            = indx_1(randperm(size(indx_1,1)));
    indx_2            = indx_2(randperm(size(indx_2,1)));

    labels_perm       = labels;
    labels_perm(indx_1(1:(length(indx_1)/2))) = 2;
    labels_perm(indx_2(1:(length(indx_2)/2))) = 1;
    cfg.design        = labels_perm;
%    

    if doqsub
    jobid{t,rep}      = qsubfeval('ft_timelockstatistics',cfg,datatmp,'memreq',1024^3,'timreq',100,'batchid',strcat('naive_',int2str(t),'_',int2str(rep)));
    else
    end
  %collect results for each time slice to avoid memory problems
  for rep2 = 20:50
        argout          = qsubget(jobid{t,rep2},'timeout',360);
        acc(t,rep2)     = argout.statistic.accuracy;
        if cfg.max_smp > 0
            trainacc(t,rep2) = argout.trainacc.statistic.accuracy;
        end
  end
end
toc



