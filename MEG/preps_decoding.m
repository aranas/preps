
%% Load data and specify variables
load('/project/3011210.01/MEG/p5_VAvsNA.mat')
timestep              = 0.1;
t                     = [timestep:timestep:2.5];
repeat                = 20;

%% prepare data
cfg                   = [];
cfg.resamplefs        = 300;
data_res              = ft_resampledata(cfg,data);

%decompose matrix and only keep 60 components
cfg                   = [];
cfg.demean            = 'yes';
cfg.scale             =  0;
cfg.method            = 'pca';
cfg.numcomponent      = 60;
data_comp             = ft_componentanalysis(cfg, data_res);

cfg                   = [];
cfg.keeptrials        = 'yes';
avg_data              = ft_timelockanalysis(cfg,data_comp);


cfg                   = [];
cfg.method            = 'crossvalidate';
cfg.mva               = 'preps_naivebayes';%{dml.standardizer dml.svm}; % 'preps_naivebayes';%;
cfg.statistic         = {'accuracy'};
cfg.type              = 'nfold'; %'bloo' only with evenly distributed classes
cfg.nfolds            = 20;
cfg.resample          = 1;% default: false; resampling for 'loo' retuns zero
cfg.design            = labels;
cfg.poolsigma         = 0;
%cfg.max_smp           = 0;
cfg.numfeat           = 250; 
%cfg.randomseed        = 5;

%loop over time & repeats
jobid                 = cell(length(t),repeat);
acc                   = zeros(length(t),repeat);
tic
for  i = 1:length(t)
      cfgt            = [];
      cfgt.latency    = [t(i)-timestep t(i)];
      datatmp         = ft_selectdata(cfgt,avg_data);
      
  for rep = 20:50
    %class             = ft_timelockstatistics(cfg,datatmp);
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
    jobid{i,rep}      = qsubfeval('ft_timelockstatistics',cfg,datatmp,'memreq',1024^3,'timreq',100,'batchid',strcat('naive_',int2str(i),'_',int2str(rep)));
   end
  %collect results for each time slice to avoid memory problems
  for rep2 = 20:50
        argout          = qsubget(jobid{i,rep2},'timeout',360);
        acc(i,rep2)     = argout.statistic.accuracy;
        if cfg.max_smp > 0
            trainacc(i,rep2) = argout.trainacc.statistic.accuracy;
        end
  end
end
toc



