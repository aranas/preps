%preps_w2v_RSA
%loads stimuli with w2v info
%selects only nouns/verbs from middle of sentence
%computes cross-correlation matrix between all (content)word embeddings
%computes cross-correlation matrix between evoked activity of all words
%computes cross-correlation matrix between evoked activity of  all final
%nouns
%computes significance of correlation between w2v similarity matrix and
%neural similarity matrix at the time words were presented and at the time
%of presenting final word of sentence.
clear all

%flags
do_source = true;
maincfg.subj = 'sub-010';
maincfg.datasuffix  = '_lp01';
maincfg.seltrig     = [112,115,122,125,212,215,222,225,113,123,213,223,118,128,218,228];

%% Load/select data & embeddings
%get filenames & trigger
root_dir    = '/project/3011210.01/MEG';
channelfile     = fullfile(root_dir,sprintf('%s_dataclean%s.mat',maincfg.subj, maincfg.datasuffix));
lcmvfile        = fullfile(root_dir,strcat(maincfg.subj,'_preps_lcmv_parc.mat'));
[seltrig, pos] = preps_help_collecttrig(maincfg.subj, maincfg.seltrig);

%load channel level data
load(channelfile,'data')
sel = ones(length(data.trialinfo),1);
sel             = sel & ismember(data.trialinfo(:,1),seltrig);
cfg             = [];
cfg.trials      = sel;
datasel         = ft_selectdata(cfg,data);
clear data

%convert to source data
if do_source
maincfg.parcel_indx = [];
num_comp = 1;
load(lcmvfile);
source_parc.filterlabel = filterlabel;
datasel = preps_sensor2parcel(datasel,source_parc,num_comp,maincfg.parcel_indx);
clear source_parc source filterlabel
end

%select word embeddings corresponding to trial
load preps_stimuli
embd = [];
pos = [];
words = {};
wrd2trl = cell(length(datasel.trialinfo),1);
for w = 1:length(datasel.trialinfo)
    wid = datasel.trialinfo(w,2);
    wnum = num2str(datasel.trialinfo(w,1));
    wnum = str2double(wnum(3));
    if ~isempty(stimuli(wid).words(wnum).w2v)
        if ~any(ismember(words,stimuli(wid).words(wnum).lemma))
            words = [words stimuli(wid).words(wnum).lemma];
            embd = [embd; stimuli(wid).words(wnum).w2v];
            pos = [pos wnum];
        end
        id = find(ismember(words,stimuli(wid).words(wnum).lemma));
        wrd2trl{id} = [wrd2trl{id} w];
    end
end
clear stimuli
nword = length(words);
[nparc, nt] = size(datasel.trial{1});

% Average over repetition of word
brain_perword = zeros(nword,nparc,nt);
for w = 1:length(words)
    trlsel = datasel.trial(wrd2trl{w});
    brain_perword(w,:,:) = mean(cat(3,trlsel{:}),3);
end

%% Create similarity matrices for data & model and Correlate

%sort for part of speech
[n,sort_id] = sort(pos);
embd = embd(sort_id,:);
brain_perword = brain_perword(sort_id,:,:);

%model RDM based on w2v
w2v_RDM = squareform(pdist(zscore(embd')','euclidean'));
figure;imagesc(w2v_RDM)
model.rdm{1} = w2v_RDM;
model.name{1} = 'allwords_embeddings';

%compute permutations
rng(5)
num_perm = 100;
perm_idx = zeros(num_perm,nword);
for r = 1:num_perm
    perm_idx(r,:) = randperm(nword);
end

%neural RDM per parcel (searchlight across time)
toverlap = 0.5;
twidth = 0.1;
begtim = datasel.time{1}(1);
endtim = datasel.time{1}(end);
endtim = endtim - twidth;
time = linspace(begtim, endtim, round(abs(begtim-endtim) ./ ...
    (twidth - toverlap * twidth)) + 1);

jobid = {};
for t = 1:length(time)
    begtim = nearest(datasel.time{1},time(t));
    endtim = nearest(datasel.time{1},time(t)+twidth);
    data = permute(brain_perword(:,:,begtim:endtim),[1,3,2]);
    
    jobid{t} = qsubfeval('preps_rsa_correlate',data,'euclidean',model,'spearman',perm_idx,1,...
        'memreq',(1024^3)*5,'timreq',60*60,'batchid',sprintf('preps_rsa_t%i_%s',t,maincfg.subj));
end
rho = zeros(length(time),nparc);
rhorand = zeros(length(time),nparc,num_perm);
for t = 2:length(time)
    try
        out = qsubget(jobid{t});
    catch
    load([jobid{t} '_output.mat'])
    out = argout{1};
    end
    rho(t,:) = out.rho;
    rhorand(t,:,:) = out.rhorand;
end

save(fullfile(root_dir,'RSA','w2v',sprintf('%s_contentwords_rho.mat',maincfg.subj)),'rho','rhorand')
%% Collect results
root_dir    = '/project/3011210.01/MEG';
load atlas_subparc374_8k

mapped_rho = zeros(10,size(atlas.parcellationlabel,1),length(time));
mapped_rhorand = zeros(10,size(atlas.parcellationlabel,1),length(time),num_perm);
for nsub = 1:10
    maincfg.subj = sprintf('sub-%.3d',nsub);%map to full source space
    load(fullfile(root_dir,'RSA','w2v',sprintf('%s_contentwords_rho.mat',maincfg.subj)))
    
    pindx = 1:length(atlas.parcellationlabel);
    pindx([1 2 188 189]) = []; %ignore medial wall parcels
    for p = 1:nparc
        indx = pindx(p);
        mapped_rho(nsub,indx,:) = rho(:,p)';
        mapped_rhorand(nsub,indx,:,:) = permute(rhorand(:,p,:),[2,1,3]);
    end
end
all_rho = squeeze(mean(mapped_rho));
all_rhorand = squeeze(mean(mapped_rhorand));
%% Stats

cfg                  = [];
cfg.connectivity     = parcellation2connmat(atlas);
cfg.tail             = 1;
cfg.clustertail      = 1;
cfg.clusterthreshold = 'nonparametric_individual';
cfg.clusteralpha     = 0.05;
cfg.feedback         = 'text';
cfg.clusterstatistic = 'maxsum';

cfg.dim = size(all_rhorand(:,:,1));
cfg.numrandomization = size(all_rhorand,3);

all_rho(all_rho==0) = nan;
all_rhorand(all_rhorand==0) = nan;
statobs  = reshape(all_rho,[],1);
statrand = reshape(all_rhorand,[],size(all_rhorand,3));

stats = clusterstat(cfg, statrand, statobs);
fn = fieldnames(stats);
for k = 1:numel(fn)
    try
        stats.(fn{k}) = reshape(stats.(fn{k}),cfg.dim);
    end
end


%% Visualize
load(fullfile('/project/3011210.01/anatomy',maincfg.subj,strcat(maincfg.subj,'_sourcemodel.mat')));
source                = [];
source.brainordinate  = atlas;
source.brainordinate.pos = sourcemodel.pos;
source.label          = atlas.parcellationlabel;
source.time           = time;
source.dimord         = 'chan_time';
source.pow            = all_rho;%.*(stats.prob<0.05);

%for exploration
cfgp                  = [];
cfgp.funparameter     = 'pow';
%cfgp.maskparameter    = source.mask;
ft_sourcemovie(cfgp, source);


% create the 'upsampling matrix', to map parcels onto the vertices of the
% cortical sheet
% x = zeros(0,1);
% y = zeros(0,1);
% for k = 1:numel(atlas.parcellationlabel)
%     sel = atlas.parcellation==k;
%     x = cat(1,x(:),find(sel));
%     y = cat(1,y(:),ones(sum(sel),1)*k);
% end
% P = sparse(x,y,ones(numel(x),1),size(atlas.pos,1),numel(atlas.parcellationlabel));
% 
% figure;ft_plot_mesh(atlas,'vertexcolor',P*mapped_rho');
% lighting gouraud; material dull;
% view([-90 0])
% camlight;
