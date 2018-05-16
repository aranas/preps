addpath(genpath('/home/language/sopara/Matlab/fieldtrip/external/dmlt'))
addpath /home/common/matlab/fieldtrip/qsub/
addpath('/home/language/sopara/Prepositionalphrases/preps/MEG')

load('/project/3011210.01/MEG/p5_alltrials_w2v')

% cfg             = [];
% cfg.resamplefs  = 300;
% data            = ft_resampledata(cfg,data);

cfg                   = [];
cfg.demean            = 'yes';
cfg.scale             =  0;
cfg.method            = 'pca';
cfg.numcomponent      = 60;
data             = ft_componentanalysis(cfg, data);

cfg             = [];
cfg.keeptrials  = 'yes';
cfg.vartrllength = 2;
avg_data        = ft_timelockanalysis(cfg,data);

%select pos
ind             = [];
for i = 1:length(allwords)
    if any(strcmp({'NN','VVFIN'},pos(i)))
        ind     = [ind i];
    end
end

feat            = feat(ind,:);
allwords        = allwords(ind);
pos             = pos(ind);
cfg             = [];
cfg.trials      = ind;
avg_data        = ft_selectdata(cfg,avg_data);

cfg             = [];
cfg.latency     = [];
cfg.method      = 'crossvalidate';
cfg.mva         = 'ridgeregression_sa';%{dml.standardizer dml.svm}; % 'preps_naivebayes';%;
cfg.statistic   = {'eval_correlation'};
cfg.type        = 'nfold'; %'bloo' only with evenly distributed classes
cfg.nfolds      = size(feat,1)/6;
cfg.design      = feat;
cfg.randomseed  = 5;
%cfg.numlambdas  = 10;
cfg.lambda      = 7.9207e+06;
cfg.vocab       = allwords;
%cfg.max_smp     = size(feat,1)-6 ;
%cfg.checksize   = 4000000; %allow to keep cv field in cfg (although large)

timestep        = 0.1;
t               = [(avg_data.time(1)+timestep):timestep:0.5];

jobid    = cell(length(t),1);
acc             = zeros(length(t),1);

for  i = 1:length(t)
    cfgt             = [];
    cfgt.latency     = [t(i)-timestep t(i)];
    datatmp          = ft_selectdata(cfgt,avg_data);
%     cfg.avgovertime = 'yes';
%     cfg.latency     = [t(i)-timestep t(i)];
    
    for rep = 1:50%:size(grid_samples,2)
        %cfg.max_smp = grid_samples(nsmp);
        
                  feat_perm     = feat(randperm(size(feat,1)),:);
                  %feat_perm     = randn(size(feat))/4;
                  cfg.design    = feat_perm;
        jobid_perm{i,rep}    = qsubfeval('ft_timelockstatistics',cfg,datatmp,'memreq',(1024^3)*20,'timreq',600,'batchid',strcat('ridge_',int2str(i),'_rep',int2str(rep),'_leave6_perm'));
    end
end

for i = 1:size(jobid_perm,1)
  for  j = 1:size(jobid_perm,2)
        argout       = load(strcat(jobid_perm{i,j},'_output.mat'));
        acc_perm(i,j)     = argout.argout{1}.statistic{1}; 
        %trainacc(i,j)= argout.argout{1}.trainacc.statistic{1};
  end
end

%inspect labels
correct_labels = class.statistic{3};
incorrect_labels = class.statistic{4};pos
count = 0;
for i = 1:9
    for j = 1:length(correct_labels)
        count = count + all(strcmp(pos_cmb{i},correct_labels{j}));
    end
    match(i) = count;
    count = 0;
end

for i = 1:9
    for j = 1:length(incorrect_labels)
        count = count + all(strcmp(pos_cmb{i},incorrect_labels{j}));
    end
    nomatch(i) = count;
    count = 0;
end

