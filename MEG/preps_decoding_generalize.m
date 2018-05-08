
function preps_decoding_generalize(subject,path,vs, test_trigger, timestep, trainingwindow, permute,repetition)

load(strcat(path,subject,'_',vs,'_generalize_2hp.mat'))

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
t        = timestep:timestep:data.time{9}(end);

%remember order of trials for appenddata
avg_data.trialinfo(:,2) = [1:length(avg_data.trialinfo)]';

cfgt             = [];
cfgt.trials       = ~ismember(avg_data.trialinfo(:,1),test_trigger);
cfgt.latency     = trainingwindow;
data_train       = ft_selectdata(cfgt,avg_data);

%class           = ft_timelockstatistics(cfg,datatmp);
% %
for rep = 1:repetition
    for j = 1:length(t) % for each time slice in test window
        cfgt              = [];
        cfgt.trials       = ismember(avg_data.trialinfo(:,1),test_trigger);
        cfgt.latency      = [t(j)-0.1 t(j)];
        data_combi        = ft_selectdata(cfgt,avg_data);
        data_combi.time   = data_train.time;
        
        cfgt              = [];
        datatmp           = ft_appenddata(cfgt,data_train,data_combi);
        %recover original trial order
        [~,idx]           = sort(datatmp.trialinfo(:,2));
        datatmp.trial     = datatmp.trial(idx);
        datatmp.trialinfo = datatmp.trialinfo(idx,:);
        
        if permute
            indx_1          = find(labels==1);
            indx_2          = find(labels==2);
            indx_1          = indx_1(randperm(size(indx_1,1)));
            indx_2          = indx_2(randperm(size(indx_2,1)));
            
            labels_perm     = labels;
            labels_perm(indx_1(1:(length(indx_1)/2))) = 2;
            labels_perm(indx_2(1:(length(indx_2)/2))) = 1;
            cfg.design      = labels_perm;
            add2name        = '_perm';
        else
            cfg.design      = labels;
            add2name        = '';
        end
        %%
        cfg.testfolds     = {find(ismember(datatmp.trialinfo(:,1),test_trigger))};
        class             = ft_timelockstatistics(cfg,datatmp);
        acc(rep,j)          = class.statistic.accuracy;
        acctrain(rep,j)     = class.trainacc.statistic.accuracy;
    end
end
prev = datatmp.cfg.previous;
argins = {subject,path,vs, test_trigger, timestep, trainingwindow, permute,repetition};
save(strcat(path,'Classification/',subject,'/','generalize2hp _',vs,'_',num2str(trainingwindow(1)),add2name),'acc','acctrain','cfg','prev','argins','-v7.3')


