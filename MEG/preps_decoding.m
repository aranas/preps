addpath(genpath('/home/language/sopara/Matlab/fieldtrip/external/dmlt'))
addpath('/home/language/sopara/Prepositionalphrases/preps/MEG')

load('/project/3011210.01/MEG/pilot_VAvsNAvsFill.mat')


cfg             = [];
cfg.keeptrials  = 'yes';
avg_data        = ft_timelockanalysis(cfg,data);


cfg             = [];
cfg.method      = 'crossvalidate_sa';
cfg.mva         = {dml.standardizer dml.naive}; %default
cfg.statistic   = {'accuracy'};
cfg.nfolds      = 20;
cfg.resample    = true;% default: false;
cfg.design      = labels;
%cfg.checksize   = 4000000; %allow to keep cv field in cfg (although large)

%loop over time
acc             = [];
for i = 0.1:0.1:1.5
cfg.latency     = [i-0.1 i];
class           = ft_timelockstatistics(cfg,avg_data);
acc             = [acc class.statistic.accuracy];
acc
end

cfg             = [];
cfg.method      = 'crossvalidate_sa';
cfg.mva         = 'ridgeregression_sa';
cfg.statistic   = {'accuracy','binomial'};
cfg.nfolds      = 5;
cfg.numlambdas  = 10;
cfg.resample    = true;% default: false;
cfg.latency     = [0.1 0.4];
cfg.design      = labels;
%cfg.checksize   = inf; %allow to keep cv field in cfg (although large)
class           = ft_timelockstatistics(cfg,avg_data)


