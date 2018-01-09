addpath(genpath('/home/language/sopara/Matlab/fieldtrip/external/dmlt'))
addpath('/home/language/sopara/Prepositionalphrases/preps/MEG')
addpath('/home/language/sopara/Prepositionalphrases/Scripts_tom/GNB')

load('/project/3011210.01/MEG/pilot_DetvsNoun.mat')


cfg             = [];
cfg.keeptrials  = 'yes';
avg_data        = ft_timelockanalysis(cfg,data);


cfg             = [];
cfg.method      = 'crossvalidate';
cfg.mva         = {dml.standardizer dml.naive}; % 'preps_naivebayes'%;
cfg.statistic   = {'accuracy'};
cfg.type        = 'nfold';
cfg.nfolds      = 5;
cfg.resample    = true;% default: false;
cfg.latency     = [0.1 0.4];
cfg.design      = labels;
%cfg.checksize   = 4000000; %allow to keep cv field in cfg (although large)

%loop over time
acc             = [];
for i = -0.5:0.1:1
cfg.latency     = [i-0.1 i];
class           = ft_timelockstatistics(cfg,avg_data);
acc             = [acc class.statistic.accuracy];
acc
end


