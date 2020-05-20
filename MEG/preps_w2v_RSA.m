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
maincfg.subj = 'sub-002';
maincfg.datasuffix  = '_lp01';
maincfg.seltrig     = [114,124,214,224,115,125,215,225,113,123,213,223];

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
wskip = [];
embd = zeros(length(datasel.trialinfo),300);
all_pos = [];
for w = 1:length(datasel.trialinfo)
    wid = datasel.trialinfo(w,2);
    wnum = num2str(datasel.trialinfo(w,1));
    wnum = str2double(wnum(3));
    
    all_pos = [all_pos wnum];
    if ~isempty(stimuli(wid).words(wnum).w2v)
        embd(w,:) = stimuli(wid).words(wnum).w2v;    
    else
        wskip = [wskip w];
    end
end
clear stimuli
%delete those trials (words) from data, for which we don't have word
%embeddings
cfg = [];
cfg.trials = true(1,length(datasel.trialinfo));
cfg.trials(wskip) = false;
datasel = ft_selectdata(cfg,datasel);
embd(wskip,:) = [];
all_pos(wskip) = [];

nword = length(datasel.trial);
[nparc, nt] = size(datasel.trial{1});

%sort for part of speech
[n,sort_id] = sort(all_pos);
embd = embd(sort_id,:);

%model RDM based on w2v
w2v_RDM = squareform(pdist(zscore(embd')','euclidean'));
figure;imagesc(w2v_RDM)

%neural RDM per parcel (across time)
data = reshape(cell2mat(datasel.trial),[nparc, nt, nword]);
data = data(:,:,sort_id);

brain_RDM = zeros(nt,nword,nword);
for p = 1:nparc
    tmp = zscore(squeeze(data(p,:,:)))';
    brain_RDM(p,:,:) = squareform(pdist(tmp,'euclidean'));
end
clear data tmp

%correlate
n_trials = 196;
rho = [];
for p = 1:nparc
    tmp_brain = squeeze(brain_RDM(p,401:end,401:end));
    tmp_model = w2v_RDM(1:200,1:200);
    rho(p) = corr(tmp_model(tril(true(n_trials,n_trials),-1)),tmp_brain(tril(true(n_trials,n_trials),-1)),'type','spearman','rows', 'complete');
end

%Visualize
%map to full source space
load atlas_subparc374_8k
pindx = 1:length(atlas.parcellationlabel);
pindx([1 2 188 189]) = []; %ignore medial wall parcels
mapped_rsm = zeros(size(atlas.parcellationlabel));
for p = 1:nparc
    indx = pindx(p);
    mapped_rho(indx) = rho(p);
end

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

figure;ft_plot_mesh(atlas,'vertexcolor',P*mapped_rho');
lighting gouraud; material dull;
view([-90 0])
camlight;
