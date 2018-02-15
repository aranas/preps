addpath(genpath('/home/language/sopara/Matlab/fieldtrip/external/dmlt'))
addpath('/home/language/sopara/Matlab/dss2_1-0')
addpath /home/common/matlab/fieldtrip/qsub/
addpath('/home/language/sopara/Prepositionalphrases/preps/MEG')
addpath('/home/language/sopara/Prepositionalphrases/Scripts_tom/GNB')
addpath(genpath('/home/language/sopara/Matlab/distributionPlot'))

load('/project/3011210.01/MEG/p2_NounvsAdj.mat')

cfg             = [];
cfg.resamplefs  = 300;
data_res            = ft_resampledata(cfg,data);

% %emulate random data
% for l = 1:size(data.trial,2)
%     data_rand{l} = randn(270,450)*5^-13;
% end
% 
% data_res.trial = data_rand;

%decompose matrix and only keep 60 components
params                = [];
params.time           = data_res.time;

cfg                   = [];
cfg.demean            = 'yes';
cfg.scale             =  0;
cfg.method            = 'pca';
cfg.numcomponent      = 60;
% cfg.cellmode          = 'yes';
% cfg.dss.algorithm     = 'pca';
% cfg.dss.denf.function = 'denoise_avg2';
% cfg.dss.denf.params   = params;
data_comp              = ft_componentanalysis(cfg, data_res);

% cfg              = [];
% data             = ft_rejectcomponent(cfg,data);

% Transform data to normality
% for i = 1:size(data_comp.trial,2)
%     for j = 1:size(data_comp.trial{1},1)
%       rank                = tiedrank(data_comp.trial{i}(j,:));
%       p                   = rank/(length(rank) + 1); %# +1 to avoid Inf for the max point
%       data_comp.trial{i}(j,:)  = norminv(p,0,1);
%     end
% end

cfg             = [];
cfg.keeptrials  = 'yes';
avg_data        = ft_timelockanalysis(cfg,data_comp);


cfg             = [];
cfg.method      = 'crossvalidate';
cfg.mva         = 'preps_naivebayes';%{dml.standardizer dml.svm}; % 'preps_naivebayes';%;
cfg.statistic   = {'accuracy'};
cfg.type        = 'bloo';
%cfg.nfolds      = 20;
cfg.resample    = 1;% default: false; resampling for 'loo' retuns zero
%vector for testfolds
cfg.design      = labels;
cfg.poolsigma   = 0;
cfg.trainacc    = 1;
%cfg.numfeat     = 300; 
%cfg.randomseed  = 5;
%cfg.checksize   = 4000000; %allow to keep cv field in cfg (although large)

%loop over time & repeats
timestep = 0.1;
numsmp   = 4:4:size(avg_data.trial,1);
t        = [timestep:timestep:1];
jobid    = cell(length(t),size(numsmp,2));

indx_1   = find(labels==1);
indx_2   = find(labels==2);



acctest      = zeros(length(t),size(numsmp,2));
acctrain      = zeros(length(t),size(numsmp,2));
tic
for  i = 1:length(t)
      cfgt             = [];
      cfgt.latency     = [t(i)-0.1 t(i)];
      datatmp        = ft_selectdata(cfgt,avg_data);
      
  for rep = 1:size(numsmp,2)
     ind              = [indx_1(randperm(size(indx_1,1),numsmp(rep)/2)); indx_2(randperm(size(indx_2,1),min(numsmp(rep)/2,size(indx_2,1))))];
     trialind         = false(1,size(avg_data.trial,1));
     trialind(ind)    = 1;
     cfgt             = [];
     cfgt.trials      = trialind;
     datatmp2          = ft_selectdata(cfgt,datatmp);
     cfg.design       = labels(trialind);
     %class            = ft_timelockstatistics(cfg,datatmp2);
      %permute labels
%     indx_1          = find(labels==1);
%     indx_2          = find(labels==2);
%     indx_1          = indx_1(randperm(size(indx_1,1)));
%     indx_2          = indx_2(randperm(size(indx_2,1)));
% 
%     labels_perm     = labels;
%     labels_perm(indx_1(1:(length(indx_1)/2))) = 2;
%     labels_perm(indx_2(1:(length(indx_2)/2))) = 1;
%    cfg.design      = labels_perm;
%     
    jobid{i,rep}    = qsubfeval('ft_timelockstatistics',cfg,datatmp2,'memreq',1024^3,'timreq',100,'batchid',strcat('naive_',int2str(i),'_',int2str(rep)));
   end
  %collect results for each time slice to avoid memory problems
  for rep2 = 1:size(numsmp,2)
        argout          = qsubget(jobid{i,rep2},'timeout',260);
        acctest(i,rep2)      = argout.statistic.accuracy;
        acctrain(i,rep2)      = argout.trainacc.statistic.accuracy;
  end
end
toc



