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
do_source = false;
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

%neural RDM per time point
data = reshape(cell2mat(datasel.trial),[nparc, nt, nword]);
data = data(:,:,sort_id);

brain_RDM = zeros(nt,nword,nword);
for t = 1:nt
    tmp = zscore(squeeze(data(:,t,:)))';
    brain_RDM(t,:,:) = squareform(pdist(tmp,'euclidean'));
end
t = nearest(datasel.time{1},0.15);
figure;imagesc(squeeze(brain_RDM(t,:,:)))

