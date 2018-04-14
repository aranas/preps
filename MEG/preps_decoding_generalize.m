
function preps_decoding_generalize(subject,path,vs, test_trigger, timestep, trainingwindow, permute,repetition)

load(strcat(path,subject,'_',vs,'_generalize.mat'))

cfg                   = [];
cfg.resamplefs        = 300;
data_res              = ft_resampledata(cfg,data);

cfg                   = [];
cfg.demean            = 'yes';
cfg.scale             =  0;
cfg.method            = 'pca';
cfg.numcomponent      = 60;
data_comp             = ft_componentanalysis(cfg, data);

cfg                   = [];
cfg.keeptrials        = 'yes';
cfg.vartrllength      = 2 ;
avg_data              = ft_timelockanalysis(cfg,data_comp);

cfg             = [];
cfg.method      = 'crossvalidate';
cfg.mva         = 'preps_naivebayes';%{dml.standardizer dml.svm}; % 'preps_naivebayes';%;
cfg.statistic   = {'accuracy'};
cfg.type        = 'split';
cfg.resample    = 1;% default: false; resampling for 'loo' retuns zero
cfg.poolsigma   = 0;
cfg.max_smp     = sum(~ismember(avg_data.trialinfo,test_trigger));
cfg.numfeat     = 500;
%cfg.randomseed  = 5;

%loop over time & repeats
t        = timestep:timestep:data.time{1}(end);
t(end+1) = data.time{1}(end);

cfgt             = [];
cfgt.trials       = ~ismember(avg_data.trialinfo,test_trigger);
cfgt.latency     = trainingwindow;
data_train       = ft_selectdata(cfgt,avg_data);

%class           = ft_timelockstatistics(cfg,datatmp);
% %
for rep = 1:repetition
    for j = 1:length(t) % for each time slice in test window
        cfgt              = [];
        cfgt.trials       = ismember(avg_data.trialinfo,test_trigger);
        cfgt.latency      = [t(j)-0.1 t(j)];
        data_combi        = ft_selectdata(cfgt,avg_data);
        data_combi.time   = data_train.time;
        
        cfgt              = [];
        datatmp           = ft_appenddata(cfgt,data_train,data_combi);
        
        for trial = 1:size(datatmp.trial,2)
            if any(ismember([112,122,212,222],datatmp.trialinfo(trial)))
                labels_train(trial)   = 1;
            elseif any(ismember([114,124,214,224],datatmp.trialinfo(trial)))
                labels_train(trial)   = 2;
            elseif any(ismember([218,228],datatmp.trialinfo(trial)))
                labels_train(trial)   = 1;
            elseif any(ismember([118,128],datatmp.trialinfo(trial)))
                labels_train(trial)   = 2;
            end
        end
        
        if permute
            indx_1          = find(labels_train==1);
            indx_2          = find(labels_train==2);
            indx_1          = indx_1(randperm(size(indx_1',1)));
            indx_2          = indx_2(randperm(size(indx_2',1)));
            
            labels_perm     = labels_train;
            labels_perm(indx_1(1:(length(indx_1)/2))) = 2;
            labels_perm(indx_2(1:(length(indx_2)/2))) = 1;
            cfg.design      = labels_perm;
            add2name        = 'perm_';
        else
            cfg.design      = labels_train;
            add2name        = '';
        end
        %%
        cfg.testfolds     = {find(ismember(datatmp.trialinfo,test_trigger))};
        class             = ft_timelockstatistics(cfg,datatmp);
        acc(rep,j)          = class.statistic.accuracy;
        acctrain(rep,j)     = class.trainacc.statistic.accuracy;
    end
end
prev = datatmp.cfg.previous;
argins = {subject,path,vs, test_trigger, timestep, trainingwindow, permute,repetition};
save(strcat(path,'Classification/',subject,'/','generalize_',vs,'_',add2name),'acc','acctrain','cfg','prev','argins','-v7.3')


