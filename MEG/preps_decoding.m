%% Binary classifier on some MEG data file and corresponding labels
% Classifiers available are:
% preps_naivebayes
% ridgeregression_sa
% {dml.standardizer dml.svm};
% {dml.standardizer dml.blogreg};

%% Set options
load preps_stimuli
root_dir     = '/project/3011210.01/MEG/';
save_dir     = '/project/3011210.01/MEG/Classification';

%Specify choices in decoding pipeline
if ~isfield(maincfg,'mode'),            maincfg.mode  = '';                end
if ~isfield(maincfg,'classifier'),      maincfg.classifier  = '';          end

if ~isfield(maincfg,'save_full'),       maincfg.save_full = false;         end
%parameters
if ~isfield(maincfg,'subj'),            maincfg.subj = 'pilot-005';        end %should be 'all' for mscca
if ~exist('suffix','var'),          suffix       = '';         end %might change later in code depending on selected options
if ~isfield(maincfg,'datasuffix'),      maincfg.datasuffix   = '';         end %can be empty for 2hp or '0.1hp'
if ~isfield(maincfg,'twidth'),          maincfg.twidth       = 0.1;        end %width of sliding time window in s
if ~isfield(maincfg,'toverlap'),        maincfg.toverlap     = 0.8;        end %how much sliding time windows will overlap in %
if ~isfield(maincfg,'time'),            maincfg.time         = 'all';      end %total time to cover
%how to treat data
if ~isfield(maincfg,'seltrig'),         maincfg.seltrig      = '';         end %default selects all
if ~isfield(maincfg,'dattype'),         maincfg.dattype      = 'sensor';   end %which data to load for decoding. Can be sensor (default), lcmv or simulate
if ~isfield(maincfg,'dopca'),           maincfg.dopca        = false;      end
if ~isfield(maincfg,'clean_muscle'),    maincfg.clean_muscle = false;      end
if ~isfield(maincfg,'resample'),        maincfg.resample     = false;      end
%how to label data
if ~isfield(maincfg,'dopretest100'),    maincfg.dopretest100 = false;                      end
if ~isfield(maincfg,'doposttest'),      maincfg.doposttest   = false;                      end
if ~isfield(maincfg,'matchlength'),     maincfg.matchlength  = false;                      end
if ~isfield(maincfg,'dow2v'),           maincfg.dow2v        = false;                      end % if false (default) dependent variable will be part of speech label
if ~isfield(maincfg,'dow2vcateg'),      maincfg.dow2vcateg   = false;                      end
if ~isfield(maincfg,'w2vdim'),          maincfg.w2vdim       = [];                         end %empty vector = do not reduce dimensionality, if number, then do reduce to specified number of dimensions
if ~isfield(maincfg,'ncluster'),        maincfg.ncluster     = 2;                          end
%how to treat time
if ~isfield(maincfg,'time_avg'),        maincfg.time_avg     = false;                      end
if maincfg.time_avg, suffix              = [suffix '_avg'];                            end
if ~isfield(maincfg,'time_concat'),     maincfg.time_concat  = false;                      end
if maincfg.time_concat, suffix              = [suffix '_concat'];                      end
%preproc trigger infofile
if ~exist('pos','var')
    [seltrig, pos] = preps_help_collecttrig(maincfg.subj, maincfg.seltrig);
else
    seltrig = maincfg.seltrig;
end
%classifier-specific parameters
switch maincfg.classifier
    
    case {'preps_naivebayes','blogreg','svm','naive'}
        %needs classes - derive from seltrig
        if ~isfield(maincfg,'statfun'), maincfg.statfun = {'accuracy','confusion'};     end
        if ~isfield(maincfg,'numfeat'), maincfg.numfeat = 250;                          end
        
    case 'ridgeregression_sa'
        if ~isfield(maincfg,'statfun'), maincfg.statfun = {'eval_compare_correlation'};         end% could be eval_compare_rank or eval_correlation
        if ~isfield(maincfg,'lambda'),  maincfg.lambda = [];                            end
        if ~isfield(maincfg,'lambdaeval'), maincfg.lambdaeval = 'mse';                     end
    case 'lda'
        if ~isfield(maincfg,'statfun'), maincfg.statfun = 'accuracy';     end
    otherwise
        warning('no known classifier is specified')
        return;
end

%mode-specific parameters
switch maincfg.mode
    case 'normal'
        if ~isfield(maincfg,'repeats'), maincfg.repeats = 50;                          end
        
    case 'general'%generalising over time
         if ~isfield(maincfg,'folds'), maincfg.folds= 20;end
         if ~isfield(maincfg,'repeats'), maincfg.repeats = 50; end
                %specify training and testing time window/trials when generalising
                %over time
         if ~isfield(maincfg,'trainwindow'), maincfg.trainwindow  = [0.3 0.4];              end
         if ~isfield(maincfg,'testtrig'), error('labels for test samples need to be specified'); end 
         seltrig    = [seltrig;maincfg.testtrig'];
         if size(maincfg.testpos,1)<size(maincfg.testpos,2), maincfg.testpos = maincfg.testpos';                end
         pos        = vertcat(pos,maincfg.testpos);
    case 'lc'
        if ~isfield(maincfg,'repeats'), maincfg.repeats = 50;                          end
    case 'tuda'
        if ~isfield(maincfg,'repeats'), maincfg.repeats = 50;                          end
    case 'mvpatoolbox'
        if ~isfield(maincfg,'repeats'), maincfg.repeats = 50;                          end 
    otherwise
        warning('no mode specified - nothing is computed')
        return;
end

%filenames
lcmvfile        = fullfile(root_dir,strcat(maincfg.subj,'_preps_lcmv_parc.mat'));
channelfile     = fullfile(root_dir,sprintf('%s_dataclean%s.mat',maincfg.subj, maincfg.datasuffix));
artfctfile      = fullfile(root_dir,strcat(maincfg.subj,'_muscle'));
subjn           = str2num(maincfg.subj(end-2:end));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% Load & Select data (Design matrix)%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch maincfg.dattype
    case {'sensor', 'lcmv'} %if using channel level data
        fprintf('retrieving channel time courses\n')
        load(channelfile,'data')
    case 'mscca'
        if ~isfield(maincfg,'mscca_concat'), maincfg.maincfg.mscca_concat = 1; end
        suffix      = [suffix '_mscca' int2str(maincfg.mscca_concat)];
        root_dir    = '/project/3011210.01/MEG/mscca';
        files       = dir(fullfile(root_dir,'mscca_parcel*.mat'));
        if maincfg.mscca_concat == 1 % concat subjects as features/predictors
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
        elseif maincfg.mscca_concat == 2%trials from all subjects are concatenated, parcels are features (similar to haxby neuron paper)
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
        elseif maincfg.mscca_concat == 3
            str = unique(maincfg.pos);
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
                cfg.channel     = comp.label(contains(comp.label,maincfg.subj));
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

if maincfg.dopretest100 %only keep trials that scored high on pre-test
    fprintf('selecting trials that scored higher than 90% on pre-test\n')
    accs        = [stimuli(1:200).acc];
    sel         = sel & ismember(data.trialinfo(:,2),find(accs>=0.9));
    suffix      = [suffix '_pretest100'];
end

if maincfg.clean_muscle %reject trials with muscle artifacts
    fprintf('rejecting trials with high frequency muscle noise\n')
    load(artfctfile);
    sel         = sel & ~ismember(data.trialinfo(:,3),noisy_trials(:,end));
    suffix      = [suffix '_cleaned'];
end

% make trial selection
if ~strcmp(maincfg.dattype,'simulate')
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
end
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
    
    if maincfg.resample
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% Create dependent variable %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~strcmp(maincfg.dattype,'simulate')
    fprintf('retrieving trial labels based on part of speech\n')
    [~,loctrial]        = ismember(datasel.trialinfo(:,1),seltrig);
    [upos, ~ , labels]  = unique(pos(loctrial));
    clear loctrial
else
    upos = {''};
end

if maincfg.doposttest % relabel samples according to individual post tests;
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

if maincfg.dow2v  % replace categorical labels with w2v info
    [w2v,datasel,labels]   = preps_getw2v(datasel);
    suffix              = [suffix '_w2v'];
end

if maincfg.dow2vcateg
    [~,~,labels]   = preps_getw2v(datasel,'ncluster',maincfg.ncluster);
    suffix                 = [suffix '_w2vcateg'];
    upos                   = {unique(labels)}; 
end

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

%%%% Further subselections based on classes
%make subselection to have same distribution of number of letters across
%classes
if maincfg.matchlength
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
if ~isempty(maincfg.classifier)
    cfgcv.method            = 'crossvalidate';
    cfgcv.mva               = maincfg.classifier;
    cfgcv.statistic         = maincfg.statfun;
    cfgcv.type              = 'nfold'; %'bloo' only with evenly distributed classes
    cfgcv.design            = labels;
end
switch maincfg.classifier
    
    case 'preps_naivebayes'
        %needs classes - derive from seltrig
        if ~strcmp(maincfg.numfeat,'all'), cfgcv.numfeat = maincfg.numfeat; end
        cfgcv.poolsigma     = 0;
        if ~isfield(maincfg,'nfolds'), maincfg.nfolds = 20; end
        cfgcv.nfolds = maincfg.nfolds;
        if length(unique(sum(labels==labels'))) ~= 1
            cfgcv.resample = 1;
        else
            cfgcv.resample = 0;
        end
    case 'blogreg'
        cfgcv.mva = {dml.standardizer dml.blogreg};
        if ~isfield(maincfg,'nfolds'), maincfg.nfolds = 20; end
        cfgcv.nfolds = maincfg.nfolds;
        if length(unique(sum(labels==labels'))) ~= 1
            cfgcv.resample = 1;
        else
            cfgcv.resample = 0;
        end
    case 'svm'
        cfgcv.mva = {dml.standardizer dml.svm};
        if ~isfield(maincfg,'nfolds'), maincfg.nfolds = 20; end
        cfgcv.nfolds = maincfg.nfolds;
        if length(unique(sum(labels==labels'))) ~= 1
            cfgcv.resample = 1;
        else
            cfgcv.resample = 0;
        end
    case 'ridgeregression_sa'
        cfgcv.lambdaeval   = maincfg.lambdaeval;
        cfgcv.lambda       = maincfg.lambda;
        cfgcv.resample     = 0;
        div                = finddivisor(size(labels,1));
        div                = div(mod(div,2)==0);
        if ~isfield(maincfg,'nfolds'), maincfg.nfolds = size(labels,1)/div(nearest(div,100)); end
    case 'lda'
        if ~isfield(maincfg,'nfolds'), maincfg.nfolds = 20; end
        cfgcv.nfolds = maincfg.nfolds;
end
fprintf('using %i folds for cross-validation\n',maincfg.nfolds)
% compute time for loop%FIXME how to treat time in tuda case?
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
elseif length(maincfg.time) == 1
    time = [maincfg.time maincfg.time + maincfg.twidth];
else
    endtim = min(cellfun(@max,datasel.time));
    time = datasel.time{1};
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

if maincfg.time_concat
        ntime = diff(nearest(datasel.time,[0 maincfg.twidth]))+1;
        labels  = repmat(labels, ntime,1);
end

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
%% if compute accuracy
switch maincfg.mode    
    case 'normal'
        fprintf('train & test model looping over time\n')
        stat         = cell(nearest(time,endtim-maincfg.twidth/2),1);
        statshuf     = cell(nearest(time,endtim-maincfg.twidth/2),1);
        param        = cell(nearest(time,endtim-maincfg.twidth/2),maincfg.nfolds,1);
        %paramshuf    = cell(nearest(time,endtim-maincfg.twidth/2),maincfg.repeats,1);
        %remove cfg field for time & memory reasons
        if isfield(datasel,'elec')
         datasel = rmfield(datasel,{'cfg','elec','grad'});
        end
        %loop over time
        %f = waitbar(0,'looping over time slices...');

        for  t = 1:length(time)
            fprintf('timeslice %u: %d to %d ms\n',t,round(time(t)*1000),round((time(t)+maincfg.twidth)*1000))
            rng('default'); % ensure same 'random' folding behaviour for each time slice.
            %waitbar(t/nearest(time,endtim-twidth),f,sprintf('timeslice %u: %d to %d ms\n repetition %d',t,round(time(t)*1000),round((time(t)+twidth)*1000),0));
            % with time dimension
            if maincfg.time_avg
                cfgcv.latency       = [time(t) time(t)+maincfg.twidth];
                cfgcv.avgovertime   = 'yes'; 
            elseif maincfg.time_concat
                seltime             = nearest(datasel.time,[time(t) time(t)+maincfg.twidth]);
               
                [n,ncomp,ntime]     = size(datasel.trial(:,:,seltime(1):seltime(2)));
                datatmp.trial       = reshape(permute(datasel.trial(:,:,seltime(1):seltime(2)) ,[1 3 2]),[n*ntime ncomp]);
                datatmp.time        = datasel.time(1);
                
                cfgcv.design        = labels;
                %labels_perm         = cellfun(@(x) repmat(x, ntime,1),labels_perm,'UniformOutput',false); 
            else
                cfgcv.latency           = [time(t) time(t)+maincfg.twidth];
            end
            out                = ft_timelockstatistics(cfgcv,datatmp);
            stat{t}            = out.statistic;
            if maincfg.save_full
                nout                    = size(out.out,1);
                param(t,1:nout)     = out.out;
            end
            %loop over repetitions
            for rep = 1:maincfg.repeats
                %waitbar(t/nearest(time,endtim-twidth)+(0.01*rep),f,sprintf('timeslice %u: %d to %d ms\n repetition %d',t,round(time(t)*1000),round((time(t)+twidth)*1000),rep));
                cfgtmp              = cfgcv;
                cfgtmp.design       = labels_perm{rep};
                outshuf             = ft_timelockstatistics(cfgtmp,datatmp);
                
                statshuf{t,rep}    = outshuf.statistic;
                
            end
        end
        %close(f)
        
        %%save results to file including timesteps
        cfgcv.time      = time;
        cfgcv.twidth    = maincfg.twidth;
        cfgcv.pca       = maincfg.dopca;
        cfgcv.vocab     = upos;
        cfgcv.trialinfo = datasel.trialinfo;
        cfgcv.previous  = datatmp.cfg;
        cfgcv.previous  = rmfield(cfgcv.previous,'previous');
        if maincfg.save_full
            cfgcv.param = param;
            %cfgcv.paramshuf = paramshuf;
            cfgcv.pca_unmixing = datasel.unmixing;
        end

        switch maincfg.classifier
            case 'preps_naivebayes'
                if strcmp(maincfg.numfeat,'all'), maincfg.numfeat = length(out.out{1}.Mu);end
                cfgcv.numfeat   = maincfg.numfeat;
                if maincfg.dow2vcateg
                    filename = fullfile(save_dir, maincfg.subj, maincfg.dattype, sprintf('nbayes_%s%s_%dfolds_%dfeats_%dcluster%s',maincfg.subj,maincfg.datasuffix,maincfg.nfolds,numfeat,length(upos{:}),suffix));
                elseif strcmp(maincfg.dattype,'lcmv') && ~isempty(maincfg.parcel_indx)
                    filename = fullfile(save_dir, maincfg.subj, maincfg.dattype, 'searchlight',sprintf('nbayes_%s%s_%dfolds_%dfeats_%s_parcel%03d%s',maincfg.subj,maincfg.datasuffix,maincfg.nfolds,maincfg.numfeat,horzcat(upos{:}),maincfg.parcel_indx,suffix));
                else
                filename = fullfile(save_dir, maincfg.subj, maincfg.dattype, sprintf('nbayes_%s%s_%dfolds_%dfeats_%s%s',maincfg.subj,maincfg.datasuffix,maincfg.nfolds,maincfg.numfeat,horzcat(upos{:}),suffix));
                end
                
            case 'ridgeregression_sa'
                filename = fullfile(save_dir, maincfg.subj, maincfg.dattype, sprintf('regress_%s%s_%dfolds_lambda%d_%s%s',maincfg.subj,maincfg.datasuffix,maincfg.nfolds,maincfg.lambda,horzcat(upos{:}),suffix));
            case 'svm'
                if strcmp(maincfg.numfeat,'all'), maincfg.numfeat = length(out.time)*length(out.label);end
                cfgcv.numfeat   = maincfg.numfeat;
                filename = fullfile(save_dir, maincfg.subj, maincfg.dattype, sprintf('svm_%s%s_%dfolds_%dfeats_%s%s',maincfg.subj,maincfg.datasuffix,maincfg.nfolds,maincfg.numfeat,horzcat(upos{:}),suffix));
            case 'blogreg'
                if strcmp(maincfg.numfeat,'all'), maincfg.numfeat = length(out.time)*length(out.label);end
                cfgcv.numfeat   = maincfg.numfeat;
                filename = fullfile(save_dir, maincfg.subj, maincfg.dattype, sprintf('blogreg_%s%s_%dfolds_%dfeats_%s%s',maincfg.subj,maincfg.datasuffix,maincfg.nfolds,cfgcv.numfeat,horzcat(upos{:}),suffix));

        end
        save(filename, 'stat','statshuf','cfgcv','-v7.3');
  
    case 'general'
        
        suffix = [suffix '_general'];
        
        stat         = cell(nearest(time,endtim-maincfg.twidth/2),1);
        statshuf     = cell(nearest(time,endtim-maincfg.twidth/2),1);
        
        %remove cfg field for time & memory reasons
        if isfield(datasel,'elec')
         datasel = rmfield(datasel,{'cfg','elec','grad'});
        end
        
        %select data to be trained on
        cfgt            = [];
        cfgt.trials     = ~ismember(datasel.trialinfo(:,1),maincfg.testtrig);
        cfgt.latency    = maincfg.trainwindow;
        data_train      = ft_selectdata(cfgt,datasel);
        
        %loop over time
        %f = waitbar(0,'looping over time slices...');
        for  t = 1:length(time)
            fprintf('timeslice %u: %d to %d ms\n',t,round(time(t)*1000),round((time(t)+maincfg.twidth)*1000))
            rng('default'); % ensure same 'random' folding behaviour for each time slice.
            %waitbar(t/nearest(time,endtim-twidth),f,sprintf('timeslice %u: %d to %d ms\n repetition %d',t,round(time(t)*1000),round((time(t)+twidth)*1000),0));
            % with time dimension
            
            %select data to test on
            cfgt            = [];
            cfgt.trials     = ismember(datasel.trialinfo(:,1),maincfg.testtrig);
            cfgt.latency    = [time(t) time(t)+maincfg.twidth];
            data_test      = ft_selectdata(cfgt,datasel);

            %combine data
            cfgt = [];
            datatmp = ft_appenddata(cfgt,data_train,data_test);
            %recover original trial order
            [~,idx] = sort(datatmp.trialinfo(:,end));
            datatmp.trial = datatmp.trial(idx);
            datatmp.trialinfo = datatmp.trialinfo(idx,:);
            datatmp.time = repmat({datatmp.time{1}},[1,size(datatmp.time,2)]);
            %set parameters for crossvalidate
            cfgcv.type = 'split';
            cfgcv.max_smp = sum(~ismember(datatmp.trialinfo(:,1),maincfg.testtrig));
            cfgcv.testfolds = {find(ismember(datatmp.trialinfo(:,1),maincfg.testtrig))};

            out                = ft_timelockstatistics(cfgcv,datatmp);
            stat{t}            = out.statistic;
            %loop over repetitions
            for rep = 1:maincfg.repeats
                %waitbar(t/nearest(time,endtim-twidth)+(0.01*rep),f,sprintf('timeslice %u: %d to %d ms\n repetition %d',t,round(time(t)*1000),round((time(t)+twidth)*1000),rep));
                cfgtmp              = cfgcv;
                cfgtmp.design       = labels_perm{rep};
                outshuf             = ft_timelockstatistics(cfgtmp,datatmp);               
                statshuf{t,rep}    = outshuf.statistic;
                
            end
        end
         %%save results to file including timesteps
        cfgcv.time      = time;
        cfgcv.twidth    = maincfg.twidth;
        cfgcv.trainwind = maincfg.trainwindow;
        cfgcv.testtrig  = maincfg.testtrig;
        cfgcv.testpos   = maincfg.testpos;
        cfgcv.pca       = maincfg.dopca;
        cfgcv.vocab     = upos;
        cfgcv.trialinfo = datasel.trialinfo;
        cfgcv.previous  = datatmp.cfg;
        cfgcv.previous  = rmfield(cfgcv.previous,'previous');
        
        if strcmp(maincfg.numfeat,'all'), maincfg.numfeat = length(out.out{1}.Mu);end
        cfgcv.numfeat   = maincfg.numfeat;
        
        trainpos = unique(maincfg.testpos);
        
        filename = fullfile(save_dir, maincfg.subj, maincfg.dattype, ...
            sprintf('nbayes_%s%s_%dfolds_%dfeats_%s%i-%i_%s',...
            maincfg.subj,maincfg.datasuffix,maincfg.nfolds,...
            maincfg.numfeat,horzcat(trainpos{:}),round(maincfg.trainwindow(1)*1000),round(maincfg.trainwindow(2)*1000),suffix));
        save(filename, 'stat','statshuf','cfgcv','-v7.3');
        
    
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


