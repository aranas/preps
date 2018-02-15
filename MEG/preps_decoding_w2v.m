addpath(genpath('/home/language/sopara/Matlab/fieldtrip/external/dmlt'))
addpath /home/common/matlab/fieldtrip/qsub/
addpath('/home/language/sopara/Prepositionalphrases/preps/MEG')

load('/project/3011210.01/MEG/p2_alltrials_w2v')

cfg             = [];
cfg.keeptrials  = 'yes';
avg_data        = ft_timelockanalysis(cfg,data);

%reduce data
feat = feat(1:800,:);
avg_data = ft_selectdata(avg_data,'rpt',[1:800]);

cfg             = [];
cfg.method      = 'crossvalidate';
cfg.mva         = 'ridgeregression_sa';%{dml.standardizer dml.svm}; % 'preps_naivebayes';%;
cfg.statistic   = {'eval_euclideandistance'};
cfg.type        = 'nfold'; %'bloo' only with evenly distributed classes
cfg.nfolds      = length(feat)/8;
cfg.design      = feat;
cfg.randomseed  = 5;
cfg.numlambdas  = 3;
%cfg.checksize   = 4000000; %allow to keep cv field in cfg (although large)

timestep = 0.01;
t        = [timestep:timestep:0.5];  

jobid    = cell(length(t),1);
acc      = zeros(length(t),1);

tic
for  i = 1:length(t)
      cfgt             = [];
      cfgt.latency     = [t(i)-timestep t(i)];
      datatmp        = ft_selectdata(cfgt,avg_data);
  
       jobid{i}    = qsubfeval('ft_timelockstatistics',cfg,datatmp,'memreq',(1024^3)*20,'timreq',6000,'batchid',strcat('ridge_',int2str(i)));
end

for  i = 1:length(t)
        argout          = qsubget(jobid{i});
        acc(i,rep2)      = argout.statistic.accuracy;
end