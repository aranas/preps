%% This script runs L2 ridge regression predicting fastText word embeddings and generalizing over time 
%% to find reactivation of previous verb/noun semantics

load preps_stimuli
root_dir     = '/project/3011210.01/MEG/';
save_dir     = '/project/3011210.01/MEG/Classification';

maincfg = [];
maincfg.mode = 'general';
maincfg.classifier = 'ridgeregression_sa';
maincfg.dow2v = 1;
maincfg.w2vdim = 5;
maincfg.statfun = {'eval_correlation'};
maincfg.lambda = 1;
maincfg.toverlap = 0.2;
maincfg.time = [0 2];
maincfg.repeats = 50;
maincfg.seltrig = {'NN', 'VVFIN', 'ADJA'};
maincfg.testtrig = {'VA'};

suffix = 'general_eval2';

%parameters
if ~isfield(maincfg,'subj'),            maincfg.subj = 'pilot-005';        end %should be 'all' for mscca
if ~isfield(maincfg,'datasuffix'),      maincfg.datasuffix   = '';         end %can be empty for 2hp or '0.1hp'
if ~isfield(maincfg,'twidth'),          maincfg.twidth       = 0.1;        end %width of sliding time window in s
if ~isfield(maincfg,'toverlap'),        maincfg.toverlap     = 0.8;        end %how much sliding time windows will overlap in %
if ~isfield(maincfg,'time'),            maincfg.time         = 'all';      end %total time to cover
%how to treat data
if ~isfield(maincfg,'seltrig'),         maincfg.seltrig      = '';         end %default selects all
if ~isfield(maincfg,'dattype'),         maincfg.dattype      = 'sensor';   end %which data to load for decoding. Can be sensor (default), lcmv or simulate
if ~isfield(maincfg,'dopca'),           maincfg.dopca        = false;      end
if ~isfield(maincfg,'clean_muscle'),    maincfg.clean_muscle = false;      end
if ~isfield(maincfg,'avgrpt'),          maincfg.avgrpt       = false;      end
if ~isfield(maincfg,'resample'),        maincfg.resample     = false;      end

if ~isfield(maincfg,'w2vdim'),          maincfg.w2vdim       = [];                         end %empty vector = do not reduce dimensionality, if number, then do reduce to specified number of dimensions

%classifier-specific parameters
if ~isfield(maincfg,'statfun'), maincfg.statfun = {'eval_compare_correlation'};         end% could be eval_compare_rank or eval_correlation
if ~isfield(maincfg,'lambda'),  maincfg.lambda = [];                            end
if ~isfield(maincfg,'lambdaeval'), maincfg.lambdaeval = 'mse';                     end

if ~isfield(maincfg,'repeats'), maincfg.repeats = 50; end

%filenames
lcmvfile        = fullfile(root_dir,strcat(maincfg.subj,'_preps_lcmv_parc.mat'));
channelfile     = fullfile(root_dir,sprintf('%s_dataclean%s.mat',maincfg.subj, maincfg.datasuffix));
artfctfile      = fullfile(root_dir,strcat(maincfg.subj,'_muscle'));
subjn           = str2num(maincfg.subj(end-2:end));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% Load & Select data (Design matrix)%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
fprintf('retrieving channel time courses\n')
load(channelfile,'data')

[seltrig, postrain] = preps_help_collecttrig(maincfg.subj, maincfg.seltrig);
[testtrig, postest] = preps_help_collecttrig(maincfg.subj, maincfg.testtrig);

sel = ones(length(data.trialinfo),1);

% make trial selection
sel             = sel & ismember(data.trialinfo(:,1),seltrig);
cfg             = [];
cfg.trials      = sel;
datasel         = ft_selectdata(cfg,data);
%fix:why does this step also slightly shift time axis?
tn = size(datasel.time{1},2);
for i = 1:length(datasel.trial)
    datasel.time{i}(1:tn) = datasel.time{1}(1:tn);
end
clear sel tn i
%clear some memory
clear data

if strcmp(maincfg.dattype,'lcmv')
    if ~isfield(maincfg,'parcel_indx'), maincfg.parcel_indx = [];       end
    if ~isempty(maincfg.parcel_indx), num_comp = 5; else, num_comp = 1; end
    
    fprintf('retrieving source time courses\n')
    if ~exist(lcmvfile,'file') == 2 %if source data doesn't exist compute from scratch
        load(channelfile,'data')
        %trigger     = [30:39,110:119,120:129,210:219,220:229];
        cfg         = [];
        cfg.trials  = ismember(data.trialinfo(:,1),cat(2,seltrig));
        data        = ft_selectdata(cfg,data);
        
        [source_parc, filterlabel, source] = preps_lcmv(maincfg.subj, data);
        save(lcmvfile,'source_parc','filterlabel','source','-v7.3')
    else% otherwise simply load
        load(lcmvfile);
    end
    source_parc.filterlabel = filterlabel;
    
    if resample
        cfg                 = [];
        cfg.resamplefs      = 150;
        cfg.detrend         = 'no';
        datasel             = ft_resampledata(cfg, datasel);
        suffix              = [suffix '_150hz'];
        
    end
    
    datasel = preps_sensor2parcel(datasel,source_parc,num_comp,maincfg.parcel_indx);
    clear source_parc source filterlabel
end

if maincfg.dopca
    fprintf('decomposing data into 60 principal components\n')
    cfgtmp                   = [];
    cfgtmp.demean            = 'yes';
    cfgtmp.scale             =  0;
    cfgtmp.method            = 'pca';
    cfgtmp.numcomponent      = 60;
    datasel                  = ft_componentanalysis(cfgtmp, datasel);
end

if maincfg.avgrpt
    %doesn't work for all triggers defined
    indxreps = [];
    for i = 1:length(datasel.trialinfo)
        id = datasel.trialinfo(i,2);
        ipos = num2str(datasel.trialinfo(i,1));
        ipos = str2num(ipos(end));
        if ~ismember(i,indxreps)
            if stimuli(id).condition == 1
                paired = find(ismember([stimuli.condition],stimuli(id).condition) & ismember([stimuli.pair_num],stimuli(id).pair_num));
                if ipos == 3 % verb repeats across sentence pairs
                    paired = find(ismember([stimuli.condition],stimuli(id).condition) & ismember([stimuli.verb_num],stimuli(id).verb_num));
                end
            elseif stimuli(id).condition == 2
                if ipos == 2; ipos = 5; elseif ipos == 5; ipos = 2; end
                paired = find(ismember([stimuli.condition],stimuli(id).condition) & ismember([stimuli.pair_num],stimuli(id).pair_num));
            end
            if stimuli(id).condition ~=3
                paired = paired(~ismember([stimuli(paired).id],id));
                paired = find(datasel.trialinfo(:,2)==paired);
                pairpos = num2str(datasel.trialinfo(paired,1));
                indx = [];
                for l = 1:size(pairpos,1)
                    if str2double(pairpos(l,3))==ipos
                        indx = paired(l);
                    end
                end
                indxreps = [indxreps indx];
                
                datasel.trial{i} = mean(cat(3,datasel.trial{indx},datasel.trial{i}),3);
            end
        end
    end  
    sel = ones(length(datasel.trialinfo),1);
    sel(indxreps) = 0;
    sel = boolean(sel);
    cfg = [];
    cfg.trials = sel;
    datasel = ft_selectdata(cfg,datasel);
    if ~dow2v
    cfg.trials = ~ismember(datasel.trialinfo(:,1),[112 122 212 222]);
    datasel     = ft_selectdata(cfg,datasel);
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% Create dependent variable %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[w2v, datasel, labels]   = preps_getw2v(datasel);
suffix              = [suffix '_w2v'];

if ~isempty(maincfg.w2vdim)
    % compute data cross-covariance matrix
    labels = bsxfun(@minus,labels,mean(labels)); %mean-center columns
    
    C = (labels'*labels)./(size(labels,1)-1);
    
    % eigenvalue decomposition (EVD)
    [V,D] = eig(C);
    
    % sort eigenvectors in descending order of eigenvalues
    d = cat(2,(1:1:size(labels,2))',diag(D));
    d = sortrows(d, -2);
    
    % return the desired number of principal components  
    unmixing = V(:,d(1:maincfg.w2vdim,1))';
    labels = unmixing*labels';
    labels = labels';
    clear C D V d unmixing
   
end

if isequal(testtrig,[119 129]') %if trying to decode verb semantics from post sentence onset.
    id = datasel.trialinfo(ismember(datasel.trialinfo(:,1),[119 129]),2);
    for i = 1:length(id)
        idx = find(ismember(datasel.trialinfo(:,2),id(i)) & ismember(datasel.trialinfo(:,1),[113,123,213,223]));
        labels_test(i,:) = labels(idx,:);
    end
    labels =  vertcat(labels,labels_test);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Set Configuration %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cfgcv = [];
cfgcv.method            = 'crossvalidate';
cfgcv.mva               = maincfg.classifier;
cfgcv.statistic         = maincfg.statfun;
cfgcv.type              = 'split';
cfgcv.nfolds            = 1;
cfgcv.design            = labels;
cfgcv.lambdaeval   = maincfg.lambdaeval;
cfgcv.lambda       = maincfg.lambda;
cfgcv.resample     = 0;

if strcmp(maincfg.time,'all')
    begtim  = min(cellfun(@min,datasel.time));
    endtim  = min(cellfun(@max,datasel.time));
    endtim = endtim - (maincfg.toverlap*maincfg.twidth);
time = linspace(begtim, endtim, round(abs(begtim-endtim) ./ ...
    (maincfg.twidth - maincfg.toverlap * maincfg.twidth)) + 1);
elseif length(maincfg.time) == 2
    begtim = maincfg.time(1);
    endtim = maincfg.time(2);
    endtim = endtim - (maincfg.toverlap*maincfg.twidth);
time = linspace(begtim, endtim, round(abs(begtim-endtim) ./ ...
    (maincfg.twidth - maincfg.toverlap * maincfg.twidth)) + 1);
else
    endtim = min(cellfun(@max,datasel.time));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% Loop over time & repeats %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% ft_timelockstatistics %%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(maincfg.dattype,'sensor')
    cfg                   = [];
    cfg.keeptrials        = 'yes';
    cfg.vartrllength      = 1;
    datasel               = ft_timelockanalysis(cfg,datasel);
end

datatmp = datasel;
%permute labels
fprintf('precomputing randomisations...\n')
labels_perm = cell(maincfg.repeats,1);

rng(5);
for rep = 1:maincfg.repeats
    if length(unique(labels))==2% if two classes do stratified permutation
        %FIXME: adapt for multiclass case.
        indx_1            = find(labels==1);
        indx_2            = find(labels==2);
        indx_1            = indx_1(randperm(size(indx_1,1)));
        indx_2            = indx_2(randperm(size(indx_2,1)));
        
        labels_perm{rep}       = labels;
        labels_perm{rep}(indx_1(1:(floor(length(indx_1)/2)))) = 2;
        labels_perm{rep}(indx_2(1:(floor(length(indx_2)/2)))) = 1;
    else
        labels_perm{rep}       = labels(randperm(size(labels,1)),:);
    end
end

%select time window to be trained on
timetrain = time(nearest(time,0.2):nearest(time,0.5));
for  t = 1:nearest(timetrain,endtim-maincfg.twidth/2)
    fprintf('timeslice %u: %d to %d ms\n',t,round(timetrain(t)*1000),round((timetrain(t)+maincfg.twidth)*1000))
    rng('default'); % ensure same 'random' folding behaviour for each time slice.
    
    cfgt             = [];
    cfgt.trials       = ismember(datasel.trialinfo(:,1),seltrig);
    cfgt.latency     = [timetrain(t) timetrain(t)+maincfg.twidth];
    data_train       = ft_selectdata(cfgt,datasel);

    for  t2 = 1:nearest(time,endtim-maincfg.twidth/2)
        fprintf('timeslice %u: %d to %d ms\n',t2,round(time(t2)*1000),round((time(t2)+maincfg.twidth)*1000))
        
        cfgt              = [];
        cfgt.trials       = ismember(datasel.trialinfo(:,1),testtrig);
        cfgt.latency      = [time(t2) time(t2)+maincfg.twidth];
        data_test        = ft_selectdata(cfgt,datasel);
        %data_test.time   = data_train.time;
        
        cfgt              = [];
        datatmp           = ft_appenddata(cfgt,data_train,data_test);
        %recover original trial order
        %[~,idx]           = sort(datatmp.trialinfo(:,end));
        %datatmp.trial     = datatmp.trial(idx);
        %datatmp.trialinfo = datatmp.trialinfo(idx,:);
        datatmp.time      = repmat({datatmp.time{1}},[1,size(datatmp.time,2)]);
        
        %set parameters
        cfgtmp            = cfgcv;
        %cfgtmp.max_smp    = length(data_train.trial);
        cfgtmp.testfolds  = {length(data_train.trial)+1:length(datatmp.trial)};
        
        out             = ft_timelockstatistics(cfgtmp,datatmp);
        stat{t,t2}         = out.statistic;
        %acctrain{t,t2}     = out.trainacc.statistic;
        
        for rep = 1:maincfg.repeats
            cfginner            = cfgtmp;
            cfginner.design     = labels_perm{rep};
            
            outshuf             = ft_timelockstatistics(cfginner,datatmp);
            statshuf{t,t2,rep}     = outshuf.statistic;
            %accshuftrain{t,t2,rep} = outshuf.trainacc.statistic;
        end
    end
end
%%save results to file including timesteps
cfgcv.timetrain   = timetrain;
cfgcv.timetest        = time;
cfgcv.twidth      = maincfg.twidth;
cfgcv.testtrigger = testtrig;
cfgcv.trialinfo = datasel.trialinfo;
cfgcv.previous  = datatmp.cfg;

filename = fullfile(save_dir, maincfg.subj, maincfg.dattype, sprintf('regress_%s%s_%dfolds_lambda%d_%s',maincfg.subj,maincfg.datasuffix,cfgcv.nfolds,cfgcv.lambda,suffix));

save(filename, 'stat','statshuf','cfgcv','-v7.3');

