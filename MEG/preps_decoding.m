%% Binary classifier on some MEG data file and corresponding labels
% Classifiers available are:
% preps_naivebayes
% {dml.standardizer dml.svm};
% {dml.standardizer dml.blogreg};

%% Set variables
load preps_stimuli
root_dir     = '/project/3011210.01/MEG/';
save_dir     = '/project/3011210.01/MEG/Classification';

%need to be specified!!
if ~exist('mode',               'var'), mode  = '';                            end
if ~exist('classifier',         'var'), classifier  = '';                      end

if ~exist('save_full',          'var'), save_full = false;                     end
%parameters
if ~exist('subj',           'var'), subj         = 'pilot-005';                end %should be 'all' for mscca
if ~exist('suffix',         'var'), suffix       = '';                         end %might change later in code depending on selected options
if ~exist('datasuffix',     'var'), datasuffix   = '';                         end %can be empty for 2hp or '0.1hp'
if ~exist('twidth',         'var'), twidth       = 0.1;                        end %width of sliding time window in s
if ~exist('toverlap',       'var'), toverlap     = 0.8;                        end %how much sliding time windows will overlap in %
if ~exist('time',           'var'), time         = 'all';                      end %total time to cover
%how to treat data
if ~exist('seltrig',        'var'), seltrig      = '';                         end %default selects all
if ~exist('dattype',        'var'), dattype      = 'sensor';                   end %which data to load for decoding. Can be sensor (default), lcmv or simulate
if ~exist('dopca',          'var'), dopca        = false;                      end
if ~exist('clean_muscle',   'var'), clean_muscle = false;                      end
if ~exist('avgrpt',         'var'), avgrpt       = false;                      end
if ~exist('resample',       'var'), resample     = false;                      end
%how to label data
if ~exist('dopretest100',   'var'), dopretest100 = false;                      end
if ~exist('doposttest',     'var'), doposttest   = false;                      end
if ~exist('matchlength',    'var'), matchlength  = false;                      end
if ~exist('dow2v',          'var'), dow2v        = false;                      end % if false (default) dependent variable will be part of speech label
if ~exist('dow2vcateg',     'var'), dow2vcateg   = false;                      end
if ~exist('ncluster',       'var'), ncluster     = 2;                          end
%how to treat time
if ~exist('time_avg',       'var'), time_avg     = false;                      end
if time_avg, suffix              = [suffix '_avg'];                            end
if ~exist('time_concat',    'var'), time_concat  = false;                      end
if time_concat, suffix              = [suffix '_concat'];                      end
%preproc trigger infofile
if ~exist('pos','var') || isempty(seltrig)
[seltrig, pos] = preps_help_collecttrig(subj, seltrig);
end
%classifier-specific parameters
switch classifier
    
    case {'preps_naivebayes','blogreg','svm','naive'}
        %needs classes - derive from seltrig
        if ~exist('statfun',    'var'), statfun = {'accuracy','confusion'};     end
        if ~exist('numfeat',    'var'), numfeat = 250;                          end
        
    case 'ridgeregression_sa'
        if ~exist('statfun',    'var'), statfun = {'eval_correlation'};         end% could be eval_rank
        if ~exist('lambda',     'var'), lambda = [];                            end
        if ~exist('lambdaeval', 'var'), lambdaeval = 'mse';                     end
    case 'lda'
        if ~exist('statfun',    'var'), statfun = 'accuracy';     end
    otherwise
        warning('no known classifier is specified')
        return;
end

%mode-specific parameters
switch mode
    case 'normal'
        if ~exist('repeats',    'var'), repeats = 50;                          end
        
    case 'general'%generalising over time
         if ~exist('folds',  'var'), folds= 20;end
         if ~exist('repeats',    'var'), repeats = 50; end
                %specify training and testing time window/trials when generalising
                %over time
         if ~exist('trainwindow','var'), trainwindow  = [0.3 0.4];              end
         if ~exist('testtrig',   'var'), error('labels for test samples need to be specified'); end 
         seltrig    = [seltrig;testtrig'];
         if size(testpos,1)<size(testpos,2), testpos = testpos';                end
         pos        = vertcat(pos,testpos);
    case 'lc'
        if ~exist('repeats',    'var'), repeats = 50;                          end
    case 'tuda'
        if ~exist('repeats',    'var'), repeats = 50;                          end
    case 'mvpatoolbox'
        if ~exist('repeats',    'var'), repeats = 50;                          end 
    otherwise
        warning('no mode specified - nothing is computed')
        return;
end

%filenames
lcmvfile        = fullfile(root_dir,strcat(subj,'_preps_lcmv_parc.mat'));
channelfile     = fullfile(root_dir,sprintf('%s_dataclean%s.mat',subj, datasuffix));
artfctfile      = fullfile(root_dir,strcat(subj,'_muscle'));
subjn           = str2num(subj(end-2:end));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% Load & Select data (Design matrix)%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch dattype
    case {'sensor', 'lcmv'} %if using channel level data
        fprintf('retrieving channel time courses\n')
        load(channelfile,'data')
    case 'mscca'
        if ~exist('mscca_concat', 'var'), mscca_concat = 1; end
        suffix      = [suffix '_mscca' int2str(mscca_concat)];
        root_dir    = '/project/3011210.01/MEG/mscca';
        files       = dir(fullfile(root_dir,'mscca_parcel*.mat'));
        if mscca_concat == 1 % concat subjects as features/predictors
            parcellabel = {};
            for f = 1:length(files)%for each parcel
                f
                load(fullfile(root_dir,files(f).name),'comp')
                cfg             = [];
                cfg.trials      = ismember(comp.trialinfo(:,1),seltrig);
                comp            = ft_selectdata(cfg,comp);
                ntrials         = size(comp.trial,2);
                [nchan,ntime]   = size(comp.trial{1});
                parceldata(f,:,:,:) = reshape(cell2mat(comp.trial),[nchan ntime ntrials]);
                parcellabel     = vertcat(parcellabel,strrep(comp.label,'mscca001',sprintf('parcel%03d',f)));
            end
            data = comp;
            for t = 1:ntrials
                data.trial{t}   = reshape(squeeze(parceldata(:,:,:,t)),[f*nchan ntime]);
                data.label      = parcellabel;
            end
        elseif mscca_concat == 2%trials from all subjects are concatenated, parcels are features (similar to haxby neuron paper)
            load atlas_subparc374_8k
            for f = 1:length(files)%for each parcel
                load(fullfile(root_dir,files(f).name),'comp')
                cfg             = [];
                cfg.trials      = ismember(comp.trialinfo(:,1),seltrig);
                comp            = ft_selectdata(cfg,comp);
                if f == 1
                    ntrials         = size(comp.trial,2);
                    [nchan,ntime]   = size(comp.trial{1});
                    data            = comp;
                    data.trial      = cell(1,size(comp.trial,2)*nchan);
                end
                for chan = 1:nchan
                    for t = 1:ntrials
                    ind = t + (ntrials*chan) - ntrials;
                    data.trial{ind} = [data.trial{ind};comp.trial{t}(chan,:)];
                    end
                end
            end
            data.time       = repmat(data.time,1,nchan);
            data.trialinfo  = repmat(data.trialinfo,nchan,1);
            data.label      = atlas.parcellationlabel;
            data.label([1 2 188 189]) = [];
        elseif mscca_concat == 3
            str = unique(pos);
            str = strcat(str{:});
            filename = sprintf('/project/3011210.01/MEG/mscca/avgdata_%s.mat',str);
            if isfile(filename)
                load(filename)
            else
                load atlas_subparc374_8k
                label = atlas.parcellationlabel;
                label([1 2 188 189]) = [];
                for f = 1:length(files)%for each parcel
                    f
                    load(fullfile(root_dir,files(f).name),'comp');
                    cfg             = [];
                    cfg.trials      = ismember(comp.trialinfo(:,1),seltrig);
                    cfg.avgoverchan = 'yes';
                    comp            = ft_selectdata(cfg,comp);
                    if f==1
                        ntrials = size(comp.trial,2);
                        data = comp;
                        data.trial = cell(1,size(comp.trial,2));
                    end
                    for t = 1:ntrials
                        data.trial{t} = [data.trial{t}; comp.trial{t}];
                    end
                end
                data.label = label;
                save(filename,'data')
            end
        elseif mscca_concat == 4 %extract aligned data for each subject separately
            load atlas_subparc374_8k
            label = atlas.parcellationlabel;
            label([1 2 188 189]) = [];
            for f = 1:length(files)%for each parcel
                f
                load(fullfile(root_dir,files(f).name),'comp')
                cfg             = [];
                cfg.trials      = ismember(comp.trialinfo(:,1),seltrig);
                cfg.channel     = comp.label(contains(comp.label,subj));
                comp         = ft_selectdata(cfg,comp);
                if f == 1
                    data            = comp;
                    data.trial      = cell(1,size(comp.trial,2));
                end
                for t = 1:size(comp.trial,2)
                    data.trial{t} = [data.trial{t};comp.trial{t}]; 
                end
            end
            data.label = label;
        end
    case 'simulate'
        fprintf('simulating channel time courses\n')
        load(channelfile,'data')
        %dependent variable == w2v
        labels = [];
        for i = 1:length(stimuli)
            for j = 1:9
                if any(strcmp(stimuli(i).words(j).pos,{'VVFIN','NN'}))
                labels = [labels; stimuli(i).words(j).w2v];
                end
            end
        end
        
        %independent variable is linear mixture of random betas + noise
        [N,q]       = size(labels);        
        [p,ttrial]  = size(data.trial{1});
        
        % random regression coefficients
        rng(5)
        beta        = rand(q,p);
        
        datasel     = data;
        datasel     = rmfield(datasel,'sampleinfo');
        datasel.trialinfo   = 1:length(labels);
        
        datasel.trial = cell(N,1); % data
        datasel.time  = cell(N,1);
        for n = 1:N
            t = 1;
            for t=1:ttrial
                if t<nearest(data.time{1},0) %if baseline
                    Y = labels(1,randperm(size(labels,2)));
                else
                    Y = labels(n,:);
                end
                datasel.trial{n}(:,t) = Y * beta(:,:) + 0.1 * randn(1,p);  
            end
            datasel.time{n} = data.time{1};
        end

        suffix       = [suffix '_simulated']; 
    case 'hyperalign'
        %load weights
    otherwise
        warning('no datatype is specified')
        return
end

sel = ones(length(data.trialinfo),1);

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

% make trial selection
if ~strcmp(dattype,'simulate')
    sel             = sel & ismember(data.trialinfo(:,1),seltrig);
    cfg             = [];
    cfg.trials      = sel;
    datasel         = ft_selectdata(cfg,data);
    %fix:why does this step also slightly shift time axis?
    tn = size(datasel.time{1},2);
    for i = 1:length(datasel.trial)
        datasel.time{i}(1:tn) = datasel.time{1}(1:tn);
    end
end
%clear some memory
clear data

if strcmp(dattype,'lcmv')
    if ~exist('parcel_indx',    'var'), parcel_indx = [];       end
    if ~isempty(parcel_indx), num_comp = 5; else, num_comp = 1; end
    
    fprintf('retrieving source time courses\n')
    if ~exist(lcmvfile,'file') == 2 %if source data doesn't exist compute from scratch
        load(channelfile,'data')
        %trigger     = [30:39,110:119,120:129,210:219,220:229];
        cfg         = [];
        cfg.trials  = ismember(data.trialinfo(:,1),cat(2,seltrig));
        data        = ft_selectdata(cfg,data);
        
        [source_parc, filterlabel, source] = preps_lcmv(subj, data);
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
    
    datasel = preps_sensor2parcel(datasel,source_parc,num_comp,parcel_indx);
    clear source_parc source filterlabel
end

if dopca
    fprintf('decomposing data into 60 principal components\n')
    cfgtmp                   = [];
    cfgtmp.demean            = 'yes';
    cfgtmp.scale             =  0;
    cfgtmp.method            = 'pca';
    cfgtmp.numcomponent      = 60;
    datasel                  = ft_componentanalysis(cfgtmp, datasel);
end

if avgrpt
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
if ~strcmp(dattype,'simulate')
    fprintf('retrieving trial labels based on part of speech\n')
    [~,loctrial]        = ismember(datasel.trialinfo(:,1),seltrig);
    [upos, ulabel , labels]  = unique(pos(loctrial));
else
    upos = {''};
end

if doposttest % relabel samples according to individual post tests;
    fprintf('relabeling trials according to individual post test\n')
    if length(upos)>2 
        warning('relabeling only makes sense when selected classes are VA & NA')
    end
    for smp = 1:size(datasel.trial,2)
        id = datasel.trialinfo(smp,2);
        if strcmp(stimuli(id).post_test(subjn).attachment,'Nomen')
            labels(smp) = find(strcmp(upos,'NA'));
        else
            labels(smp) = find(strcmp(upos,'VA'));
        end
    end
    suffix      = [suffix '_posttest'];
end

if dow2v  % replace categorical labels with w2v info
    [w2v,datasel,labels]   = preps_getw2v(datasel);
    suffix              = [suffix '_w2v'];
end

if dow2vcateg
    [w2v,datasel,labels]   = preps_getw2v(datasel,'ncluster',ncluster);
    suffix                 = [suffix '_w2vcateg'];
    upos                   = {unique(labels)}; 
end

%%%% Further subselections based on classes
%make subselection to have same distribution of number of letters across
%classes
if matchlength
    ntrials = length(datasel.trialinfo);
    n = zeros(1,ntrials);
    for i = 1:ntrials
        tmppos = num2str(datasel.trialinfo(i,1));
        n(i) = length(stimuli(datasel.trialinfo(i,2)).words(str2double(tmppos(end))).word{1});
    end
    indx1 = find(labels==1);
    indx2 = find(labels==2);
    figure;
    h1 = histogram(n(indx1));hold on
    h2 = histogram(n(indx2));
    
    minVal = min(h1.Values,h2.Values);
    nletter = [3:11];
    rng(9)
    sel = [];
    for l = 1:length(nletter) % for each length
        indl = find(n==nletter(l));
        tmp1 = intersect(indx1,indl)';
        tmp2 = intersect(indx2,indl)';
        sel = [sel tmp1(randperm(length(tmp1),minVal(l))) tmp2(randperm(length(tmp2),minVal(l)))];
    end
    cfg         = [];
    cfg.trials  = false(1,length(datasel.trialinfo));
    cfg.trials(sel) = 1;
    datasel     = ft_selectdata(cfg,datasel);
    labels      = labels(cfg.trials);
    
    suffix      = [suffix '_matchlength'];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Set Configuration %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cfgcv = [];
if ~isempty(classifier)
    cfgcv.method            = 'crossvalidate';
    cfgcv.mva               = classifier;
    cfgcv.statistic         = statfun;
    cfgcv.type              = 'nfold'; %'bloo' only with evenly distributed classes
    cfgcv.design            = labels;
end
switch classifier
    
    case 'preps_naivebayes'
        %needs classes - derive from seltrig
        if ~strcmp(numfeat,'all'), cfgcv.numfeat = numfeat; end
        cfgcv.poolsigma     = 0;
        if ~exist('nfolds','var'), nfolds = 20; end
        cfgcv.nfolds = nfolds;
        if length(unique(sum(labels==labels'))) ~= 1
            cfgcv.resample = 1;
        else
            cfgcv.resample = 0;
        end
    case 'blogreg'
        cfgcv.mva = {dml.standardizer dml.blogreg};
        if ~exist('nfolds','var'), nfolds = 20; end
        cfgcv.nfolds = nfolds;
        if length(unique(sum(labels==labels'))) ~= 1
            cfgcv.resample = 1;
        else
            cfgcv.resample = 0;
        end
    case 'svm'
        cfgcv.mva = {dml.standardizer dml.svm};
        if ~exist('nfolds','var'), nfolds = 20; end
        cfgcv.nfolds = nfolds;
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
        if ~exist('nfolds','var'), cfgcv.nfolds = size(labels,1)/div(nearest(div,20)); end
    case 'lda'
        if ~exist('nfolds','var'), nfolds = 20; end
        cfgcv.nfolds = nfolds;
end
fprintf('using %i folds for cross-validation\n',cfgcv.nfolds)
% compute time for loop%FIXME how to treat time in tuda case?
if strcmp(time,'all')
    begtim  = min(cellfun(@min,datasel.time));
    endtim  = min(cellfun(@max,datasel.time));
    endtim = endtim - (toverlap*twidth);
time = linspace(begtim, endtim, round(abs(begtim-endtim) ./ ...
    (twidth - toverlap * twidth)) + 1);
elseif length(time) == 2
    begtim = time(1);
    endtim = time(2);
    endtim = endtim - (toverlap*twidth);
time = linspace(begtim, endtim, round(abs(begtim-endtim) ./ ...
    (twidth - toverlap * twidth)) + 1);
else
    endtim = min(cellfun(@max,datasel.time));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% Loop over time & repeats %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% ft_timelockstatistics %%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(dattype,'sensor')
    cfg                   = [];
    cfg.keeptrials        = 'yes';
    cfg.vartrllength      = 1;
    datasel               = ft_timelockanalysis(cfg,datasel);
end

datatmp = datasel;
%permute labels
fprintf('precomputing randomisations...\n')
labels_perm = cell(repeats,1);

if time_concat
        ntime = diff(nearest(datasel.time,[0 twidth]))+1;
        labels  = repmat(labels, ntime,1);
end

rng(5);
for rep = 1:repeats
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
%% if compute accuracy
switch mode    
    case 'normal'
        fprintf('train & test model looping over time\n')
        stat         = cell(nearest(time,endtim-twidth/2),1);
        statshuf     = cell(nearest(time,endtim-twidth/2),1);
        param        = cell(nearest(time,endtim-twidth/2),repeats,1);
        paramshuf    = cell(nearest(time,endtim-twidth/2),repeats,1);
        %remove cfg field for time & memory reasons
        if isfield(datasel,'elec')
         datasel = rmfield(datasel,{'cfg','elec','grad','sampleinfo'});
        end
        %loop over time
        %f = waitbar(0,'looping over time slices...');

        for  t = 1:nearest(time,endtim-twidth/2)
            fprintf('timeslice %u: %d to %d ms\n',t,round(time(t)*1000),round((time(t)+twidth)*1000))
            rng('default'); % ensure same 'random' folding behaviour for each time slice.
            %waitbar(t/nearest(time,endtim-twidth),f,sprintf('timeslice %u: %d to %d ms\n repetition %d',t,round(time(t)*1000),round((time(t)+twidth)*1000),0));
            % with time dimension
            if time_avg
                cfgcv.latency       = [time(t) time(t)+twidth];
                cfgcv.avgovertime   = 'yes'; 
            elseif time_concat
                seltime             = nearest(datasel.time,[time(t) time(t)+twidth]);
               
                [n,ncomp,ntime]     = size(datasel.trial(:,:,seltime(1):seltime(2)));
                datatmp.trial       = reshape(permute(datasel.trial(:,:,seltime(1):seltime(2)) ,[1 3 2]),[n*ntime ncomp]);
                datatmp.time        = datasel.time(1);
                
                cfgcv.design        = labels;
                %labels_perm         = cellfun(@(x) repmat(x, ntime,1),labels_perm,'UniformOutput',false); 
            else
                cfgcv.latency           = [time(t) time(t)+twidth];
            end
            %loop over repetitions
            for rep = 1:repeats
                %waitbar(t/nearest(time,endtim-twidth)+(0.01*rep),f,sprintf('timeslice %u: %d to %d ms\n repetition %d',t,round(time(t)*1000),round((time(t)+twidth)*1000),rep));
                out                = ft_timelockstatistics(cfgcv,datatmp);
                
                cfgtmp              = cfgcv;
                cfgtmp.design       = labels_perm{rep};
                outshuf             = ft_timelockstatistics(cfgtmp,datatmp);
                
                stat{t,rep}        = out.statistic;
                statshuf{t,rep}    = outshuf.statistic;
                
                if save_full
                    nout                    = size(out.out,1);
                    param(t,rep,1:nout)     = out.out;
                    paramshuf(t,rep,1:nout) = outshuf.out;
                end
            end
        end
        %close(f)
        
        %%save results to file including timesteps
        cfgcv.time      = time;
        cfgcv.twidth    = twidth;
        cfgcv.vocab     = upos;
        cfgcv.trialinfo = datasel.trialinfo;
        cfgcv.previous  = datatmp.cfg;
        cfgcv.previous  = rmfield(cfgcv.previous,'previous');
        if save_full
            cfgcv.param = param;
            cfgcv.paramshuf = paramshuf;
        end

        switch classifier
            case 'preps_naivebayes'
                if strcmp(numfeat,'all'), numfeat = length(out.out{1}.Mu);end
                cfgcv.numfeat   = numfeat;
                if dow2vcateg
                    filename = fullfile(save_dir, subj, dattype, sprintf('nbayes_%s%s_%dfolds_%dfeats_%dcluster%s',subj,datasuffix,cfgcv.nfolds,numfeat,length(upos{:}),suffix));
                elseif strcmp(dattype,'lcmv') && ~isempty(parcel_indx)
                    filename = fullfile(save_dir, subj, dattype, 'searchlight',sprintf('nbayes_%s%s_%dfolds_%dfeats_%s_parcel%03d%s',subj,datasuffix,cfgcv.nfolds,numfeat,horzcat(upos{:}),parcel_indx,suffix));
                else
                filename = fullfile(save_dir, subj, dattype, sprintf('nbayes_%s%s_%dfolds_%dfeats_%s%s',subj,datasuffix,cfgcv.nfolds,numfeat,horzcat(upos{:}),suffix));
                end
                
            case 'ridgeregression_sa'
                filename = fullfile(save_dir, subj, dattype, sprintf('regress_%s%s_%dfolds_lambda%d_%s%s',subj,datasuffix,cfgcv.nfolds,cfgcv.lambda,horzcat(upos{:}),suffix));
            case 'svm'
                filename = fullfile(save_dir, subj, dattype, sprintf('svm_%s%s_%dfolds_%dfeats_%s%s',subj,datasuffix,cfgcv.nfolds,numfeat,horzcat(upos{:}),suffix));
            case 'blogreg'
                filename = fullfile(save_dir, subj, dattype, sprintf('blogreg_%s%s_%dfolds_%dfeats_%s%s',subj,datasuffix,cfgcv.nfolds,numfeat,horzcat(upos{:}),suffix));

        end
        save(filename, 'stat','statshuf','cfgcv','-v7.3');
        
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
        
        suffix              = [suffix '_general'];
        %select time window to be trained on
        cfgt             = [];
        cfgt.trials       = ~ismember(datasel.trialinfo(:,1),testtrig);
        cfgt.latency     = trainwindow;
        data_train       = ft_selectdata(cfgt,datasel);
        
        for t = 1:nearest(time,endtim-twidth/2)
            cfgt              = [];
            cfgt.trials       = ismember(datasel.trialinfo(:,1),testtrig);
            cfgt.latency      = [time(t) time(t)+twidth];
            data_test        = ft_selectdata(cfgt,datasel);
            %data_test.time   = data_train.time;
            
            cfgt              = [];
            datatmp           = ft_appenddata(cfgt,data_train,data_test);
            %recover original trial order
             [~,idx]           = sort(datatmp.trialinfo(:,end));
             datatmp.trial     = datatmp.trial(idx);
             datatmp.trialinfo = datatmp.trialinfo(idx,:);
             datatmp.time      = repmat({datatmp.time{1}},[1,size(datatmp.time,2)]);
             
            %set parameters
            cfgtmp            = cfgcv;
            cfgtmp.type       = 'split';
            cfgtmp.max_smp    = sum(~ismember(datatmp.trialinfo(:,1),testtrig));
            cfgtmp.testfolds  = {find(ismember(datatmp.trialinfo(:,1),testtrig))};
            
            out             = ft_timelockstatistics(cfgtmp,datatmp);
            stat(t)         = out.statistic.accuracy;
            acctrain(t)     = out.trainacc.statistic.accuracy;
            
            parfor rep = 1:repeats
                cfginner            = cfgtmp;
                cfginner.design     = labels_perm{rep};
                
                outshuf             = ft_timelockstatistics(cfginner,datatmp);
                statshuf(t,rep)     = outshuf.statistic.accuracy;
                accshuftrain(t,rep) = outshuf.trainacc.statistic.accuracy;
            end
        end
        %%save results to file including timesteps
        cfgcv.time        = time;
        cfgcv.twidth      = twidth;
        cfgcv.trainwindow = trainwindow;
        cfgcv.testtrigger = testtrig;
        cfgcv.testlabel   = testpos;
        cfgcv.vocab     = upos;
        cfgcv.trialinfo = datasel.trialinfo;
        cfgcv.previous  = datatmp.cfg;
        if strcmp(numfeat,'all'), numfeat = length(out.out{1}.Mu);end
        cfgcv.numfeat   = numfeat;
                
        switch classifier
            case 'preps_naivebayes'
                if dow2vcateg
                    filename = fullfile(save_dir, subj, dattype, sprintf('nbayes_%s%s_%dfolds_%dfeats_%dcluster%s',subj,datasuffix,cfgcv.nfolds,numfeat,length(upos{:}),suffix));
                elseif strcmp(dattype,'lcmv') && ~isempty(parcel_indx)
                    filename = fullfile(save_dir, subj, dattype, 'searchlight',sprintf('nbayes_%s%s_%dfolds_%dfeats_%s_parcel%03d%s',subj,datasuffix,cfgcv.nfolds,numfeat,horzcat(upos{:}),parcel_indx,suffix));
                else
                filename = fullfile(save_dir, subj, dattype, sprintf('nbayes_%s%s_%dfolds_%dfeats_%s%s',subj,datasuffix,cfgcv.nfolds,numfeat,horzcat(upos{:}),suffix));
                end
                
            case 'ridgeregression_sa'
                filename = fullfile(save_dir, subj, dattype, sprintf('regress_%s%s_%dfolds_lambda%d_%s%s',subj,datasuffix,cfgcv.nfolds,cfgcv.lambda,horzcat(upos{:}),suffix));
            case 'svm'
                filename = fullfile(save_dir, subj, dattype, sprintf('svm_%s%s_%dfolds_%dfeats_%s%s',subj,datasuffix,cfgcv.nfolds,numfeat,horzcat(upos{:}),suffix));
            case 'blogreg'
                filename = fullfile(save_dir, subj, dattype, sprintf('blogreg_%s%s_%dfolds_%dfeats_%s%s',subj,datasuffix,cfgcv.nfolds,numfeat,horzcat(upos{:}),suffix));

        end
        save(filename, 'stat','statshuf','cfgcv','-v7.3');
        
        %% if do Hidden markov model as described in Vidaurre et al. 2018
    case 'tuda'
        for  t = 1:nearest(time,endtim-twidth)
            fprintf('timeslice %u: %d to %d ms\n',t,round(time(t)*1000),round((time(t)+twidth)*1000))
            rng('default'); % ensure same 'random' folding behaviour for each time slice.
            cfgt = [];
            cfgt.latency           = [time(t) time(t)+twidth];
            datasel2    = ft_selectdata(cfgt,datasel);
            
            K                       = 4;
            [N p ttrial]            = size(datasel2.trial);
            datatmp                 = reshape(permute(datasel2.trial,[3 1 2]),[ttrial*N p]);
            datatmp(any(isnan(datatmp),2),:) = [];
            T                       = repmat(length(datasel2.time),[N 1]);
            
            options.K               = K;
            options.parallel_trials = 1;
            [tuda,Gamma]            = tudatrain (datatmp,labels,T,options);
            
            cfgcv.Gamma             = Gamma;
            for rep = 1:repeats
                out                     = ft_timelockstatistics(cfgcv,datasel2);
                
                
                cfgtmp              = cfgcv;
                cfgtmp.design       = labels_perm{rep};
                outshuf             = ft_timelockstatistics(cfgtmp,datasel2);
                
                stat{t,rep}        = out.statistic;
                statshuf{t,rep}    = outshuf.statistic;
                
                if ~isfield(out,'out')
                end
            end
        end
        
        suffix = [suffix '_tuda'];
        %%save results to file including timesteps
        cfgcv.time      = time;
        cfgcv.twidth    = twidth;
        cfgcv.vocab     = upos;
        cfgcv.Gamma     = Gamma;
        cfgcv.trialinfo = datasel.trialinfo;
        filename = fullfile(save_dir, subj, sprintf('classacc_%s%s_%dfolds_%dfeats_%s%s',subj,datasuffix,cfgcv.nfolds,numfeat,horzcat(upos{:}),suffix));
        save(filename, 'stat','statshuf','cfgcv');
    
    
    case 'mvpatoolbox'
        cfg = [];
        cfg.trials = true(1,length(datasel.trial));
        for i = 1:length(datasel.trial)
            if any(isnan(datasel.trial{i}(:)))
                cfg.trials(i) = false;
            end
        end
        datatmp = ft_selectdata(cfg,datasel);
        labels = labels(cfg.trials);
        
        rng('default');
        stat = cell(1,nearest(time,endtim-twidth/2));
        for  t = 1:nearest(time,endtim-twidth/2)
            cfg = [];
            cfg.latency       = [time(t) time(t)+twidth];
            
            cfg.avgovertime   = 'no';
            cfg.mvpa.cattim   = 1;
            cfg.method = 'mvpa';
            cfg.mvpa.classifier = 'svm';%goes fast for lda, but takes forever with svm
            cfg.mvpa.metric = 'accuracy';
            cfg.mvpa.k = nfolds;
            cfg.mvpa.repeat = 5;
            cfg.mvpa.stratify = 1;
            %cfg.timextime = 'yes';
            cfg.design = labels';
            stat{t} = ft_timelockstatistics(cfg,datatmp);
        end
        figure
        hold on
        for  t = 1:nearest(time,endtim-twidth/2)
            mv_plot_result(stat{t}.mvpa,stat{t}.time)
        end
end


