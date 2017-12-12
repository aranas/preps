function [single_tlck,featv,pos] = singletrial_chop(subjectname)

% load in raw subject data and selection of words/trials
% Preprocess data and select extracted trials

%read in trial selection from txt file

% file = trialselection;
% fileID =fopen(file,'r')
% formatSpec = '%d %d %s';
% stimsel = textread(fileID,formatSpec);

f   = mous_db_getfilename(subjectname, 'meg_ds_task');
if numel(f)>1
    ext = cell(1,numel(f));
    for k = 1:numel(f)
        ext{k} = sprintf('_pt%s',num2str(k));
    end
else
    ext{1} = '';
end

data = cell(1,numel(f));
mask = cell(1,numel(f));
for k = 1:numel(f)
    trl        = mous_defineTrial(f{k}, 0, 0, 'trialfun_visual_sentence');
    trign = [1 5 2 6]; %take sentence condition only
    sel = ismember(trl(:,5),trign);
    trl = trl(sel,:);
end

cfg            = [];
cfg.dataset    = f{k};
cfg.trl        = trl;
cfg.continuous = 'yes';
cfg.channel    = 'MEG';
cfg.dftfilter  = 'yes';
cfg.dftfreq    = [50 100 150 200 250 300];
cfg.padding    = 2;
cfg.demean     = 'yes';
data{k}        = ft_preprocessing(cfg);

% create a mask variable that codes for the artifacts
mous_db_getdata(subjectname, ['meg_artifact_cfg',ext{k}]);
if ~isempty(cfgjump.artfctdef.zvalue.artifact)
    % take half the data padding length for preprocessing
    cfgjump.artfctdef.zvalue.artifact(:,1) = cfgjump.artfctdef.zvalue.artifact(:,1)-1200*5;
    cfgjump.artfctdef.zvalue.artifact(:,2) = cfgjump.artfctdef.zvalue.artifact(:,2)+1200*5;
end
trlnew = mous_artifact_remove(trl, f{k}, {cfgeog1 cfgeog2 cfgjump cfgmuscle});

dum    = false(1,max(trl(:,2)));
for kk = 1:size(trlnew,1)
    dum(trlnew(kk,1):trlnew(kk,2))=true;
end

mask{k} = cell(1,numel(data{k}.trial));
for kk = 1:numel(data{k}.trial)
    mask{k}{1,kk} = dum(data{k}.sampleinfo(kk,1):data{k}.sampleinfo(kk,2));
end

% keep track of the gradiometer info
sens(k)    = data{k}.grad;
weights(k) = numel(data{k}.trial);


if numel(f)>1
    cutoffn = data{1}.trialinfo(end,1);
    data = ft_appenddata([], data{:});
    mask = cat(2, mask{:});
    data.grad = ft_average_sens(sens, 'weights', weights);
else
    data = data{1};
    mask = mask{1};
end

% Correct for projector delay in visual condition
if strcmp(subjectname(1), 'V')
    data.time = cellfun(@(x)x-0.036,data.time,'Un',0);
end

% set data to nan where artefact
for i = 1:length(data.trial)
    data.trial{i}(mask{i}) = 0;
end

trlallwords = [];
% load in word onset smpl information via trialfun_visual_words
for k = 1:numel(f)
    trlallwords{k}        = mous_defineTrial(f{k}, 0, 0, 'trialfun_visual_word');
end
if numel(f)>1
    indxkeep = ~(trlallwords{1}(:,4) > cutoffn);
    trlallwords{1} = trlallwords{1}(indxkeep',:);
    trlallwords = cat(1, trlallwords{:});
else
    trlallwords = trlallwords{1};
end

tmpindx = ismember(trlallwords(:,5),trign);
trlallwords = trlallwords(tmpindx,:);

indx = cell(1,numel(data.trial));
pre = cell(1,numel(data.trial));
pst = cell(1,numel(data.trial));
indxzero = find(trlallwords(:,6) == 0);

for k = 1:numel(data.trial)
    if k < numel(indxzero)
        idx = [indxzero(k):(indxzero(k+1)-1)];
        onsets = trlallwords(idx,6);
    else
        idx = [indxzero(k):length(trlallwords)];
        onsets = trlallwords(idx,6);
    end
    if all(isfinite(onsets))
        for m = 1:numel(onsets)
            indx{k}(m) = onsets(m)+nearest(data.time{k},0);
            pre{k}(m)  = indx{k}(m)-120; pre{k}(m) = max(pre{k}(m),1);
            pre{k}(m)  = indx{k}(m)-pre{k}(m);
            if m<numel(onsets)
                pst{k}(m) = (onsets(m+1)+nearest(data.time{k},0))-indx{k}(m);
            else
                pst{k}(m) = numel(data.time{k})-indx{k}(m);
            end
            pst{k} = min(pst{k}, 600);
        end
    end
end


for k = 1:numel(pre)
    pre{k} = pre{k}(:);
    pst{k} = pst{k}(:);
    indx{k} = indx{k}(:);
end

params.tr = indx;
params.pre = pre;
params.pst = pst;
params.mask = mask;
params.computenew = 0;
params.demean = 'prezero';
params.covariance = 1;
X.x = 1;
[~,~,tmpt]=denoise_avg_sa(params,data.trial,X);

% exclude trials with missing data
indmissing = [];
for i = 1:length(tmpt)
    if isempty(tmpt{i})
        indmissing = [indmissing i];
    elseif any(any(isnan(tmpt{i})))
        indmissing = [indmissing i];
    end
end
indmissing = unique(indmissing);
trlallwords(indmissing,:) = [];
tmpt(indmissing) = [];

single_tlck = [];
single_tlck.label = data.label;
single_tlck.dimord = 'chan_time';
single_tlck.time   = repmat({((1:721)-120)./1200},[size(tmpt),1]);
single_tlck.trials = tmpt;
single_tlck.grad = data.grad;
single_tlck.grad = trlallwords;

