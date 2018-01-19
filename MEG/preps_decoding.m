addpath(genpath('/home/language/sopara/Matlab/fieldtrip/external/dmlt'))
addpath /home/common/matlab/fieldtrip/qsub/
addpath('/home/language/sopara/Prepositionalphrases/preps/MEG')
addpath('/home/language/sopara/Prepositionalphrases/Scripts_tom/GNB')
addpath(genpath('home/language/sopara/Matlab/distributionPlot'))

load('/project/3011210.01/MEG/p2_NounvsAdj.mat')


cfg             = [];
cfg.keeptrials  = 'yes';
avg_data        = ft_timelockanalysis(cfg,data);

indx_1          = find(labels==1);
indx_2          = find(labels==2);
labels_perm     = labels;
labels_perm(indx_1(1:(length(indx1)/2))) = 2;
labels_perm(indx_2(1:(length(indx2)/2))) = 1;



cfg             = [];
cfg.method      = 'crossvalidate';
cfg.mva         = 'preps_naivebayes';%{dml.standardizer dml.naive}; % 'preps_naivebayes';%;
cfg.statistic   = {'accuracy'};
cfg.type        = 'nfold';
cfg.nfolds      = 5;
%cfg.resample    = 1;% default: false; resampling for 'loo' retuns zero
%vector for testfolds
cfg.design      = labels_perm;
cfg.poolsigma   = 0;
%cfg.checksize   = 4000000; %allow to keep cv field in cfg (although large)

%loop over time & repeats
timestep = 0.1;
t        = [timestep:timestep:0.5];
repeat   = 100;
jobid    = cell(length(t),repeat);

acc      = zeros(length(t),repeat);
tic
for i = 1:length(t)
cfgt             = [];
cfgt.latency     = [t(i)-0.1 t(i)];
datatmp        = ft_selectdata(cfgt,avg_data);

  for rep = 1:repeat
    %class           = ft_timelockstatistics(cfg,datatmp);
    jobid{i,rep}    = qsubfeval('ft_timelockstatistics',cfg,datatmp,'memreq',(1024^3)*2,'timreq',600,'batchid',strcat('naive_',int2str(i),'_',int2str(rep)));
  end
  %collect results for each time slice to avoid memory problems
  for rep2 = 1:repeat
        argout          = qsubget(jobid{i,rep2},'timeout',260);
        acc(i,rep2)      = argout.statistic.accuracy;
  end
end
toc
acc_perm_5f = acc;


