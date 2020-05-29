%preps_w2v_RSA
%loads stimuli with fastText info
%computes distance matrix between all word embeddings for selected
%words
%computes distance matrix between evoked activity of all selected words

clear all

%parameters
do_source = true;
datasuffix  = '_lp01';
%seltrig contains the selection of trigger values for both brain data and
%any amount of models to be tested, such that seltrig{1} = brain data and
%seltrig{2:n} = triggers for model comparison.
all_seltrig     = {[112,115,122,125,212,215,222,225,113,123,213,223,118,128,218,228]};
toverlap = 0.2;
twidth = 0.1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% For each subject compute & save %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% correlation between model RDM and neural RDM %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
jobid = {};
for nsub = 1:10
    subj = sprintf('sub-%.3d',nsub);%map to full source space
    
    %% Load/select data & embeddings
    
    root_dir    = '/project/3011210.01/MEG';
    channelfile     = fullfile(root_dir,sprintf('%s_dataclean%s.mat',subj, datasuffix));
    lcmvfile        = fullfile(root_dir,strcat(subj,'_preps_lcmv_parc.mat'));
    load(channelfile,'data')
    
    for i = 1:length(all_seltrig)
        
        [seltrig{i}, pos{i}] = preps_help_collecttrig(subj, all_seltrig{i});
        sel = ones(length(data.trialinfo),1);
        sel             = sel & ismember(data.trialinfo(:,1),seltrig{i});
        cfg             = [];
        cfg.trials      = sel;
        datasel{i}         = ft_selectdata(cfg,data);
        
    end
    clear data
    upos = uniqueStrCell(pos{1});
    suffix = strcat(upos{:});
    
    %convert to source data
    if do_source
        parcel_indx = [];
        num_comp = 1;
        load(lcmvfile);
        source_parc.filterlabel = filterlabel;
        for i = 1:length(datasel)
            datasel{i} = preps_sensor2parcel(datasel{i},source_parc,num_comp,parcel_indx);
        end
        clear source_parc source filterlabel
    end
    
    %select word embeddings corresponding to selected trials
    load preps_stimuli
    embd = [];
    PoS = [];
    wids = [];
    words = {};
    wrd2trl = cell(length(datasel{end}.trialinfo),1);
    for i = 1:max(length(datasel)-1,1)
        for w = 1:length(datasel{i}.trialinfo)
            wid = datasel{i}.trialinfo(w,2);
            wnum = num2str(datasel{i}.trialinfo(w,1));
            wnum = str2double(wnum(3));
            if ~isempty(stimuli(wid).words(wnum).w2v)
                if ~any(ismember(words,stimuli(wid).words(wnum).lemma))
                    words = [words stimuli(wid).words(wnum).lemma];
                    embd = [embd; stimuli(wid).words(wnum).w2v];
                    PoS = [PoS wnum];
                    wids = [wids wid];
                end
                id = find(ismember(words,stimuli(wid).words(wnum).lemma));
                wrd2trl{id} = [wrd2trl{id} w];
            end
        end
        all_pos{i} = PoS;
        all_wids{i} = wids;
        
        %model RDM based on embeddings
        model_RDM = squareform(pdist(zscore(embd')','euclidean'));
        %figure;imagesc(model_RDM)
        model.rdm{i} = model_RDM;
        upos = uniqueStrCell(pos{i});
        model.name{i} = strcat(upos{:});
    end
    clear stimuli model_RDM
    
    nword = length(words);
    [nparc, nt] = size(datasel{end}.trial{1});
    % Average over repeated words
    brain_perword = zeros(nword,nparc,nt);
    for w = 1:length(words)
        trlsel = datasel{end}.trial(wrd2trl{w});
        brain_perword(w,:,:) = mean(cat(3,trlsel{:}),3);
    end
    clear trlsel
    %% Create similarity matrices for data and Correlate with models  
    %neural RDM per parcel (searchlight across time) & correlate
    begtim = datasel{end}.time{1}(1);
    endtim = datasel{end}.time{1}(end);
    endtim = endtim - twidth;
    time = linspace(begtim, endtim, round(abs(begtim-endtim) ./ ...
        (twidth - toverlap * twidth)) + 1);
    
    for t = 1:length(time)
        begtim = nearest(datasel{end}.time{1},time(t));
        endtim = nearest(datasel{end}.time{1},time(t)+twidth);
        data = permute(brain_perword(:,:,begtim:endtim),[1,3,2]);
        
        jobid{nsub,t} = qsubfeval('preps_rsa_correlate',data,'euclidean',model,'spearman',[],1,...
            'memreq',(1024^3)*5,'timreq',60*30,'batchid',sprintf('preps_rsa_t%i_%s_%s',t,subj,suffix));
    end
    allinfo{nsub}.pos = all_pos;
    allinfo{nsub}.wid = all_wids;
    allinfo{nsub}.model = model;
end
%Collect results
for nsub = 1:10
    subj = sprintf('sub-%.3d',nsub);%map to full source space
 
    rho = [];
    for t = 1:length(time)
        load([jobid{nsub,t} '_output.mat'])
        out = argout{1};            
        for i = 1:length(out)
            rho(i,t,:) = out(i).rho;
        end
    end
    modelinfo = allinfo{nsub};
    save(fullfile(root_dir,'RSA','w2v',sprintf('%s_%s_rho.mat',subj,suffix)),'rho','modelinfo','time')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% Collect data of all subjects %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% compute significance & plot results %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Collect results
clear all
root_dir    = '/project/3011210.01/MEG';
load atlas_subparc374_8k

mapped_rho = [];
for nsub = 1:10
    subj = sprintf('sub-%.3d',nsub);%map to full source space
    load(fullfile(root_dir,'RSA','w2v',sprintf('%s_ADJANNVVFIN_rho.mat',subj)))
    
    pindx = 1:length(atlas.parcellationlabel);
    pindx([1 2 188 189]) = []; %ignore medial wall parcels
    for i = 1:size(rho,1) %for every model
        for p = 1:size(rho,3)
            indx = pindx(p);
            mapped_rho(i,nsub,indx,:) = squeeze(rho(i,:,p));
        end
    end
end
for i = 1:size(mapped_rho,1)
avg_rho{i} = squeeze(mean(squeeze(mapped_rho(1,:,:,:))));
end

mapped_rho(:,:,[1 2 188 189],:) = nan;
%% Stats
Nsub = 10;
modeli = 1;
load(fullfile('/project/3011210.01/anatomy',subj,strcat(subj,'_sourcemodel.mat')));
source                = [];
source.brainordinate  = atlas;
source.brainordinate.pos = sourcemodel.pos;
source.label          = atlas.parcellationlabel;
source.time           = time;
%source.dimord         = 'chan_time';


for nsub = 1:Nsub
    all_sources{nsub} = source;
    all_sources{nsub}.pow = squeeze(mapped_rho(modeli,nsub,:,:));
    all_bsl{nsub} = source;
    all_bsl{nsub}.pow = zeros(size(all_sources{nsub}.pow));
end

cfg = [];
cfg.method = 'montecarlo';
cfg.statistic = 'ft_statfun_depsamplesT';
cfg.parameter = 'pow';
cfg.correctm = 'cluster';
cfg.inside = true(1,length(atlas.parcellationlabel));
cfg.inside([1 2 188 189]) = false;
cfg.neighbours = 'biosemi16_neighb.mat'; %dummy neighbourhood structure
cfg.connectivity     = parcellation2connmat(atlas);
cfg.alpha = 0.05;
cfg.tail = 1;
cfg.design(1,1:2*Nsub) = [ones(1,Nsub) 2*ones(1,Nsub)];
cfg.design(2,1:2*Nsub) = [1:Nsub 1:Nsub];
cfg.numrandomization = 1000;
cfg.ivar = 1;
cfg.uvar = 2;
stat = ft_timelockstatistics(cfg,all_sources{:},all_bsl{:});


%% Visualize
%for exploration
source.pow = avg_rho{modeli}.*stat.mask; 

cfgp                  = [];
cfgp.funparameter     = 'pow';
%cfgp.maskparameter    = source.mask;
ft_sourcemovie(cfgp, source);


% create the 'upsampling matrix', to map parcels onto the vertices of the
% cortical sheet
x = zeros(0,1);
y = zeros(0,1);
for k = 1:numel(atlas.parcellationlabel)
    sel = atlas.parcellation==k;
    x = cat(1,x(:),find(sel));
    y = cat(1,y(:),ones(sum(sel),1)*k);
end
P = sparse(x,y,ones(numel(x),1),size(atlas.pos,1),numel(atlas.parcellationlabel));


%figure;ft_plot_mesh(atlas,'vertexcolor',P*squeeze(mapped_rho');
lighting gouraud; material dull;
view([-90 0])
camlight;
