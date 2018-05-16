function acc = preps_decoding_w2v_generalize(subject,path,label_trigger, test_trigger, timestep, trainingwindow, permute,repetition)

load(strcat(path,subject,'_','alltrials_w2v.mat'))

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

cfgt            = [];
cfgt.trials     = ~ismember(avg_data.trialinfo(:,1),test_trigger);
cfgt.latency    = trainingwindow;
data_train      = ft_selectdata(cfgt,avg_data);
%temporarily remove label vectors for test items to match data_train
feat            = feat(cfgt.trials,:);
allwords        = allwords(cfgt.trials);
pos             = pos(cfgt.trials);

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
data_train        = ft_selectdata(cfg,data_train);

%add back new label vectors
ind_feat         = ismember(data_train.trialinfo(:,1),label_trigger);
feat            = [feat;feat(ind_feat,:)];

cfg             = [];
cfg.latency     = [];
cfg.method      = 'crossvalidate';
cfg.mva         = 'ridgeregression_sa';
cfg.statistic   = {'eval_correlation'};
cfg.type        = 'nfold';
cfg.nfolds      = 1;
cfg.design      = feat;
cfg.lambda      = 7.9207e+06;
cfg.testsize    = 4;
%cfg.vocab       = pos;
%cfg.max_smp     = size(feat,1)-6 ;

t               = [(avg_data.time(1)+timestep):timestep:2.5];
jobid           = cell(length(t),1);
acc             = zeros(length(t),1);

for  i = 1:length(t)
    cfgt              = [];
    cfgt.trials       = ismember(avg_data.trialinfo(:,1),test_trigger);
    cfgt.latency      = [t(i)-timestep t(i)];
    data_combi        = ft_selectdata(cfgt,avg_data);
    data_combi.time   = data_train.time;
    
    cfgt              = [];
    datatmp           = ft_appenddata(cfgt,data_train,data_combi);
    for rep = 1:repetition
        if permute
            feat_perm     = feat(randperm(size(feat,1)),:);
            cfg.design    = feat_perm;
            add2name      = '_perm';
        else
            cfg.design    = feat;
            add2name      = '';
        end
        cfg.testfolds     = {find(ismember(datatmp.trialinfo(:,1),test_trigger))};
        class             = ft_timelockstatistics(cfg,datatmp);
        acc(i,rep)        = class.statistic{1};
        %jobid_perm{i,repetition}    = qsubfeval('ft_timelockstatistics',cfg,datatmp,'memreq',(1024^3)*15,'timreq',600,'batchid',strcat('ridge_',int2str(i),'_rep',int2str(rep),'_leave6'));
    end
end
prev = datatmp.cfg.previous;
argins = {subject,path,label_trigger, test_trigger, timestep, trainingwindow, permute,repetition};
save(strcat(path,'Classification/',subject,'/','generalize_w2v_',num2str(trainingwindow(1)),add2name),'acc','acctrain','cfg','prev','argins','-v7.3')


%inspect labels
% correct_labels = class.statistic{3};
% incorrect_labels = class.statistic{4};pos
% count = 0;
% for i = 1:9
%     for j = 1:length(correct_labels)
%         count = count + all(strcmp(pos_cmb{i},correct_labels{j}));
%     end
%     match(i) = count;
%     count = 0;
% end
%
% for i = 1:9
%     for j = 1:length(incorrect_labels)
%         count = count + all(strcmp(pos_cmb{i},incorrect_labels{j}));
%     end
%     nomatch(i) = count;
%     count = 0;
% end
%
