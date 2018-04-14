addpath(genpath('/home/language/sopara/Matlab/fieldtrip/external/dmlt'))
addpath('/home/language/sopara/Matlab/dss2_1-0')
addpath /home/common/matlab/fieldtrip/qsub/
addpath('/home/language/sopara/Prepositionalphrases/preps/MEG')
addpath('/home/language/sopara/Prepositionalphrases/Scripts_tom/GNB')
addpath(genpath('/home/language/sopara/Matlab/distributionPlot'))

load('/project/3011210.01/MEG/p5_VAvsNA.mat')

cfg             = [];
cfg.resamplefs  = 300;
data_res            = ft_resampledata(cfg,data);

%decompose matrix and only keep 60 components
cfg                   = [];
cfg.demean            = 'yes';
cfg.scale             =  0;
cfg.method            = 'pca';
cfg.numcomponent      = 60;
data_comp              = ft_componentanalysis(cfg, data_res);

cfg             = [];
cfg.keeptrials  = 'yes';
avg_data        = ft_timelockanalysis(cfg,data_comp);

cfg             = [];
cfg.method      = 'crossvalidate';
cfg.mva         = {'dml.standardizer' 'dml.blogreg(''scale'',0.0001)'}; % 'preps_naivebayes';%;
cfg.statistic   = {'accuracy'};
cfg.type        = 'nfold';
cfg.nfolds      = 20;
cfg.resample    = 1;% default: false; resampling for 'loo' retuns zero
%vector for testfolds
cfg.design      = labels;
%cfg.poolsigma   = 0;
%cfg.numfeat     = 500; 
%cfg.randomseed  = 5;
%cfg.checksize   = 4000000; %allow to keep cv field in cfg (although large)

%loop over time & repeats
timestep = 0.1;
t        = [timestep:timestep:2];
m        = length(labels);

grid_samples = [2:4:m-m/20];

jobid    = cell(length(t),size(grid_samples,2));

acctest      = zeros(length(t),size(grid_samples,2));
acctrain      = zeros(length(t),size(grid_samples,2));
tic
for  i = 1:length(t)
      cfgt             = [];
      cfgt.latency     = [t(i)-timestep t(i)];
      datatmp        = ft_selectdata(cfgt,avg_data);
      
  for nsmp = 1:size(grid_samples,2)
      cfg.max_smp = grid_samples(nsmp);
%       indx_1          = find(labels==1);
%       indx_2          = find(labels==2);
%       indx_1          = indx_1(randperm(size(indx_1,1)));
%       indx_2          = indx_2(randperm(size(indx_2,1)));
% 
%       labels_perm     = labels;
%       labels_perm(indx_1(1:(length(indx_1)/2))) = 2;
%       labels_perm(indx_2(1:(length(indx_2)/2))) = 1; 
%       cfg.design      = labels_perm;
    %argout = ft_timelockstatistics(cfg,datatmp);
     jobid{i,nsmp}    = qsubfeval('ft_timelockstatistics',cfg,datatmp,'memreq',1024^3,'timreq',600,'batchid',strcat('lc_',int2str(i),'_',int2str(grid_samples(nsmp)),'_samples'));
   end
  %collect results for each time slice to avoid memory problems
  for nsmp2 = 1:size(grid_samples,2)
        argout          = qsubget(jobid{i,nsmp2},'timeout',1000);
        acctest(i,nsmp2)      = argout.statistic.accuracy;
        acctrain(i,nsmp2)      = argout.trainacc.statistic.accuracy;
  end
end
toc