%% Binary classifier on some MEG data file and corresponding labels
% Classifiers available are:
% preps_naivebayes
% {dml.standardizer dml.svm};
% {dml.standardizer dml.blogreg};

%% Set variables

root_dir     = '/project/3011210.01/MEG/';
save_dir     = '/project/3011210.01/MEG/Classification';

%parameters
if ~exist('subj',           'var'), subj         = 'pilot-005';                end
if ~exist('suffix',         'var'), suffix       = '';                         end %might change later in code depending on selected options
if ~exist('datasuffix',     'var'), datasuffix   = '';                         end %can be empty for 2hp or '0.1hp'
if ~exist('twidth',         'var'), twidth       = 0.1;                        end %width of sliding time window in s
if ~exist('toverlap',       'var'), toverlap     = 0.5;                        end %how much sliding time windows will overlap in %
if ~exist('time',           'var'), time         = 'all';                      end %total time to cover
%data selection & preprocessing options
if ~exist('seltrig',        'var'), seltrig      = '';                         end %default selects all
if ~exist('dolcmv',         'var'), dolcmv       = false;                      end
if ~exist('dopca',          'var'), dopca        = false;                      end
if ~exist('clean_muscle',   'var'), clean_muscle = false;                      end
if ~exist('dopretest100',   'var'), dopretest100 = false;                      end
if ~exist('doposttest',     'var'), doposttest   = false;                      end
if ~exist('dow2v',          'var'), dow2v        = false;                      end % if false (default) dependent variable will be part of speech label

%preproc trigger info
[seltrig, pos] = preps_help_collecttrig(subj, seltrig);

%classifier-specific parameters
if ~exist('classifier',      'var'), classifier  = '';                          end
switch classifier
    
    case 'preps_naivebayes'
        %needs classes - derive from seltrig
        if ~exist('statfun',    'var'), statfun = {'accuracy'};                 end
        if ~exist('numfeat',    'var'), numfeat = 250;                          end
        
    case 'ridgeregression_sa'
        if ~exist('statfun',    'var'), statfun = {'eval_correlation'};         end% could be eval_rank
        if ~exist('lambda',     'var'), lambda = [];                            end
        if ~exist('lambdaeval', 'var'), lambdaeval = 'mse';                     end
    otherwise
        warning('no classifier is specified')
        return;
end

%mode-specific parameters
if ~exist('mode',               'var'), mode  = '';                             end
switch mode
    case 'normal'
        if ~exist('repeats',    'var'), repeats = 100;                          end
        %     case 'general'%generalising over time
        %         %FIXME not fully integrated yet!!!
        %         if ~exist('folds',  'var'), folds= 20;end
        %         %specify training and testing time window/trials when generalising
        %         %over time
        %         if ~exist('trainwindow','var'), trainwindow  = [0.3 0.4];              end
        %         if ~exist('testtrig',   'var'), testtrig     = horzcat(trigger{6:7});  end
        %         fprintf(strcat('generalising to events ',repmat(' %d',size(testtrig)),'\n'),testtrig)
        %         if ~any(ismember(testtrig,seltrig)) && ~exist('testlabel','var')
        %             error('labels for test samples need to be specified')
        %         else
        %             fprintf(strcat('labeling events as ',repmat(' %d',size(testlabel)),'\n'),testlabel)
        %         end
        
    case 'lc'
        if ~exist('repeats',    'var'), repeats = 100;                          end
    case 'tuda'
        if ~exist('repeats',    'var'), repeats = 100;                          end
    otherwise
        warning('no mode specified - nothing is computed')
        return;
end
%filenames
lcmvfile        = fullfile(root_dir,strcat(subj,'_preps_lcmv_parc.mat'));
channelfile     = fullfile(root_dir,sprintf('%s_dataclean%s.mat',subj, datasuffix));
artfctfile      = fullfile(root_dir,strcat(subj,'_muscle'));
subjn           = str2num(subj(end-2:end));

load preps_stimuli
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% Load & Select data (Design matrix)%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if dolcmv %if using source level data
    fprintf('retrieving source time courses\n')
    if exist(lcmvfile) == 2 %if precomputed source data exist, simply load
        load(lcmvfile);
        data        = source_parc;
        clear source_parc
    else%otherwise compute from scratch
        load(channelfile,'data')
        cfg         = [];
        cfg.trials  = ismember(data.trialinfo(:,1),cat(2,trigger{:}));
        data        = ft_selectdata(cfg,data);
        
        [source_parc, filterlabel, source] = preps_lcmv(subj, data);
        save(lcmvfile,'source_parc','filterlabel','source','-v7.3')
    end
else %if using channel level data
    fprintf('retrieving channel time courses\n')
    load(channelfile,'data')
end

% make trial selection
sel         = ismember(data.trialinfo(:,1),seltrig);

if dopretest100 %only keep trials that scored high on pre-test
    fprintf('selecting trials that scored higher than 90% on pre-test\n')
    accs        = [stimuli(1:200).acc];
    sel         = sel & ismember(data.trialinfo(:,2),find(accs>=0.9));
    suffix      = [suffix '_pretest100'];
end

if clean_muscle %reject trials with muscle artifacts
    fprintf('rejecting trials with high frequency muscle noise\n')
    load(artfctfile);
    sel         = sel & ~ismember(data.trialinfo(:,3),noisy_trials(:,end));
    suffix      = [suffix '_cleaned'];
end

cfg             = [];
cfg.trials      = sel;
datasel         = ft_selectdata(cfg,data);

%clear some memory
clear data

if dopca
    fprintf('decomposing data into 60 principal components\n')
    cfgtmp                   = [];
    cfgtmp.demean            = 'yes';
    cfgtmp.scale             =  0;
    cfgtmp.method            = 'pca';
    cfgtmp.numcomponent      = 60;
    datasel                  = ft_componentanalysis(cfgtmp, datasel);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% Create dependent variable %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('retrieving trial labels based on part of speech\n')
[~,loctrial]        = ismember(datasel.trialinfo(:,1),seltrig);
[upos, ulabel , labels]  = unique(pos(loctrial));

if doposttest % relabel samples according to individual post tests;
    fprintf('relabeling trials according to individual post test\n')
    if length(upos)>2 || contains({'VA','NA'},upos)
        warning('relabeling only makes sense when selected classes are VA & NA')
    end
    for smp = 1:size(datasel.trial,2)
        id = datasel.trialinfo(smp,2);
        if strcmp(stimuli(id).post_test(subjn).attachment,'Nomen')
            labels(smp) = ulabel(find(strcmp(upos,'NA')));
        else
            labels(smp) = ulabel(find(strcmp(upos,'VA')));
        end
    end
    suffix      = [suffix '_posttest'];
end

if dow2v  % replace categorical labels with w2v info
    fprintf('replacing categorical labels with continuous word embedding values'\n)
    labels = zeros(size(datasel.trial,2),300);
    indsel = false(1,size(datasel.trial,2));
    for smp = 1:size(datasel.trial,2)
        id          = datasel.trialinfo(smp,2);
        trig        = num2str(datasel.trialinfo(smp,1));
        w2vtmp      = stimuli(id).words(str2num(trig(end))).w2v;
        if size(w2vtmp,2)~=0
            labels(smp,:) = stimuli(id).words(str2num(trig(end))).w2v;
            indsel(smp) = 1;
        end
    end
    %only keep trials for which w2v is known
    labels          = labels(indsel,:);
    cfg             = [];
    cfg.trials      = indsel;
    datasel         = ft_selectdata(cfg,datasel);
    suffix          = [suffix '_w2v'];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Set Configuration %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isempty(classifier)
cfgcv                   = [];
cfgcv.method            = 'crossvalidate';
cfgcv.mva               = classifier;
cfgcv.statistic         = statfun;
cfgcv.type              = 'nfold'; %'bloo' only with evenly distributed classes
cfgcv.design            = labels;
end
switch classifier
    
    case 'preps_naivebayes'
        %needs classes - derive from seltrig
        cfgcv.numfeat       = numfeat;
        cfgcv.poolsigma     = 0;
        if ~exist('nfolds','var'), cfgcv.nfolds = 20; end
        if length(unique(sum(labels==labels'))) ~= 1
            cfgcv.resample = 1;
        else
            cfgcv.resample = 0;
        end
        
    case 'ridgeregression_sa'
        cfgcv.lambdaeval   = lambdaeval;
        cfgcv.lambda       = lambda;
        cfgcv.resample     = 0;
        div                = divisors(size(labels,1));
        div                = div(mod(div,2)==0);
        if ~exist('nfolds','var'), cfgcv.nfolds = size(labels,1)/div(6); end
end
fprintf('using %i folds for cross-validation\n',cfgcv.nfolds) 
% compute time for loop%FIXME how to treat time in tuda case?
if strcmp(time,'all')
    begtim  = min(cellfun(@min,datasel.time));
    endtim  = min(cellfun(@max,datasel.time));
else
    begtim = time(1);
    endtim = time(2);
end

time = linspace(begtim, endtim, round(abs(begtim-endtim) ./ ...
    (twidth - toverlap * twidth)) + 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% Loop over time & repeats %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% ft_timelockstatistics %%%%%%%%%%%%%%%%%%%%%%%%%

cfg                   = [];
cfg.keeptrials        = 'yes';
cfg.vartrllength      = 1;
datasel               = ft_timelockanalysis(cfg,datasel);

%% if compute accuracy
switch mode
    
    case 'normal'
        
        %% Loop over time & repeats
        acc     = zeros(length(time),repeats);
        accshuf = zeros(length(time),repeats);
        for  t = 1:(length(time)-1)
            cfgt            = [];
            cfgt.latency    = [time(t) time(t)+twidth];
            datatmp         = ft_selectdata(cfgt,datasel);
            rng('default'); % ensure same 'random' behaviour for each time slice.
            for rep = 1:repeats
                out               = ft_timelockstatistics(cfgcv,datatmp);
                switch classifier
                    case 'preps_naivebayes'
                        acc(t,rep)   = out.statistic.accuracy;
                    case 'ridgeregression_sa'
                        acc.acc(t,rep)   = out.statistic{1};
                        acc.model(t,rep,:,:) = out.model;
                        acc.out{t,rep} = out.out;
                end
                %permute labels
                if ulabel==2% if two classes do stratified permutation
                    %FIXME: adapt for multiclass case.
                    indx_1            = find(labels==1);
                    indx_2            = find(labels==2);
                    indx_1            = indx_1(randperm(size(indx_1,1)));
                    indx_2            = indx_2(randperm(size(indx_2,1)));
                    
                    labels_perm       = labels;
                    labels_perm(indx_1(1:(length(indx_1)/2))) = 2;
                    labels_perm(indx_2(1:(length(indx_2)/2))) = 1;
                else
                    labels_perm       = labels(randperm(size(labels,1)),:);
                end
                cfgtmp            = cfgcv;
                cfgtmp.design     = labels_perm;
                outshuf           = ft_timelockstatistics(cfgtmp,datatmp);
                switch classifier
                    case 'preps_naivebayes'
                        accshuf(t,rep)   = outshuf.statistic.accuracy;
                    case 'ridgeregression_sa'
                        cfg.lambda   = lambda;
                        accshuf.acc(t,rep)   = outshuf.statistic{1};
                        accshuf.model(t,rep,:,:) = out.model;
                end
            end
        end
        %%save results to file including timesteps
        cfgcv.timeinfo = time;
        cfgcv.twidth   = twidth;
        switch classifier
            case 'preps_naivebayes'
                filename = fullfile(save_dir, subj, sprintf('classacc_%s%s_%dfolds_%dfeats_%s%s',subj,datasuffix,cfgcv.nfolds,numfeat,horzcat(upos{:}),suffix));
                save(filename, 'acc','accshuf','cfgcv');
            case 'ridgeregression_sa'
                filename = fullfile(save_dir, subj, sprintf('classacc_%s%s_%dfolds_lambda%d_%s%s',subj,datasuffix,cfgcv.nfolds,cfgcv.lambda,horzcat(upos{:}),suffix));
                save(filename, 'acc','accshuf','cfgcv');
        end
        
    case 'lc'
        N         = length(labels);
        grid_smp  = [2:4:N-N/20];
        
        %% Loop over time & repeats
        acctest                   = zeros(length(tsteps),length(grid_smp));
        acctrain                  = zeros(length(tsteps),length(grid_smp));
        for t = 1:(length(time)-1)
            cfgt            = [];
            cfgt.latency    = [time(t) time(t)+twidth];
            datatmp         = ft_selectdata(cfgt,datasel);
            rng('default'); % ensure same 'random' behaviour for each time slice.
            parfor nsmp = 1:size(grid_smp,2)
                cfgtmp             = cfgcv; %need to assign variable within parfor loop
                cfgtmp.max_smp     = grid_smp(nsmp);
                out                = ft_timelockstatistics(cfgtmp,datatmp);
                acctest(t,nsmp)    = out.statistic.accuracy;
                acctrain(t,nsmp)   = out.trainacc.statistic.accuracy;
            end
        end
        %%save results to file including timesteps
        cfgcv.timeinfo = time;
        cfgcv.twidth   = twidth;
        filename = fullfile(save_dir, subj, sprintf('classlc_%s%s_%dfolds_%dfeats_%s%s',subj,datasuffix,folds,numfeat,horzcat(classes{:}),suffix));
        save(filename, 'acctest','acctrain','cfgcv');
        
        %% if do generalize over time
    case 'general'
        
        %remember order of trials for appenddata
        datasel.trialinfo(:,end+1) = [1:length(datasel.trialinfo)]';
        %select time window to be trained on
        cfgt             = [];
        cfgt.trials       = ~ismember(datasel.trialinfo(:,1),testtrig);
        cfgt.latency     = trainwindow;
        data_train       = ft_selectdata(cfgt,datasel);
        
        for t = 1:(length(time)-1)
            cfgt              = [];
            cfgt.trials       = ismember(datasel.trialinfo(:,1),testtrig);
            cfgt.latency      = [time(t) time(t)+twidth];
            data_test        = ft_selectdata(cfgt,datasel);
            data_test.time   = data_train.time;
            
            cfgt              = [];
            datatmp           = ft_appenddata(cfgt,data_train,data_test);
            %recover original trial order
            [~,idx]           = sort(datatmp.trialinfo(:,end));
            datatmp.trial     = datatmp.trial(idx);
            datatmp.trialinfo = datatmp.trialinfo(idx,:);
            %set parameters
            cfgtmp            = cfgcv;
            cfgtmp.type       = 'split';
            cfgtmp.max_smp    = sum(~ismember(datatmp.trialinfo(:,1),testtrig));
            cfgtmp.testfolds  = {find(ismember(datatmp.trialinfo(:,1),testtrig))};
            
            out             = ft_timelockstatistics(cfgtmp,datatmp);
            acc(t)      = out.statistic.accuracy;
            acctrain(t) = out.trainacc.statistic.accuracy;
            
            parfor rep = 1:repeats
                indx_1            = find(labels==1);
                indx_2            = find(labels==2);
                indx_1            = indx_1(randperm(size(indx_1,1)));
                indx_2            = indx_2(randperm(size(indx_2,1)));
                
                labels_perm       = labels;
                labels_perm(indx_1(1:(length(indx_1)/2))) = 2;
                labels_perm(indx_2(1:(length(indx_2)/2))) = 1;
                cfginner            = cfgtmp;
                cfginner.design     = labels_perm;
                
                outshuf             = ft_timelockstatistics(cfginner,datatmp);
                accshuf(t,rep)      = outshuf.statistic.accuracy;
                accshuftrain(t,rep) = outshuf.trainacc.statistic.accuracy;
            end
        end
        %%save results to file including timesteps
        cfgtmp.timeinfo    = time;
        cfgtmp.twidth      = twidth;
        cfgtmp.trainwindow = trainwindow;
        cfgtmp.testtrigger = testtrig;
        cfgtmp.testlabel   = testlabel;
        prev = datatmp.cfg.previous;
        testpos = pos(cellfun(@(x) any(ismember(x,testtrig)),trigger));
        
        filename = fullfile(save_dir, subj, sprintf('classgeneral_%s%s_%dfeats_%sto%s%s',subj,datasuffix,numfeat,horzcat(classes{:}),horzcat(testpos{:}),suffix));
        save(filename, 'acc','acctrain','cfgtmp','prev');
        
        filename = fullfile(save_dir, subj, sprintf('classgeneral_%s%s_%dfeats_%sto%s%s_shuf',subj,datasuffix,numfeat,horzcat(classes{:}),horzcat(testpos{:}),suffix));
        save(filename, 'accshuf','accshuftrain','cfgtmp','prev');
        
        %% if do Hidden markov model as described in Vidaurre et al. 2018
    case 'tuda'
        cfgcv.constant          = 0;
        K                       = 4;
        [N p ttrial]            = size(datasel.trial);
        datatmp                 = reshape(permute(datasel.trial,[3 1 2]),[ttrial*N p]);
        datatmp(any(isnan(datatmp),2),:) = [];
        T                       = repmat(length(datasel.time),[N 1]);
        
        options.K               = K;
        options.parallel_trials = 1;
        [tuda,Gamma]            = tudatrain (datatmp,labels,T,options);
        
        cfgcv.Gamma             = Gamma;
        out                     = ft_timelockstatistics(cfgcv,datasel);
        
end



