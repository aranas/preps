%This script loads the mscca-aligned signals from all subjects and divides
%each component into final noun-attached words and final verb-attached
%words. It then computes the between subjects correlation matrix and compares time-resolved
%correlations of noun-attached words and correlations of verb-attached
%words with correlations of between noun and verb attached words. 

%init
if ~exist('parcel_indx','var'), parcel_indx = 1;    end

root_dir    = '/project/3011210.01/MEG';
files       = dir(fullfile(root_dir,'mscca'));
files       = files(3:end);

%load parcel data
load(fullfile(root_dir,'mscca',files(parcel_indx).name))

%separate verb-attached and noun-attached final words
cfg = [];
cfg.trials = ismember(comp.trialinfo(:,1),[129,119]);
compNA = ft_selectdata(cfg,comp);
cfg.trials = ismember(comp.trialinfo(:,1),[229,219]);
compVA = ft_selectdata(cfg,comp);

%sort words according to sentence ID for both attachments separately
[~, Ni] = sort(compNA.trialinfo(:,2));
[~, Vi] = sort(compVA.trialinfo(:,2));
compNA.trial = compNA.trial(Ni);
compNA.trialinfo = compNA.trialinfo(Ni,:);
compVA.trial = compVA.trial(Vi);
compVA.trialinfo = compVA.trialinfo(Vi,:);

%add attachment info to labels and concatenate
compNA.label = strcat(compNA.label,'NA');
compVA.label = strcat(compVA.label,'VA');
cfg = [];
comp = ft_appenddata(cfg,compVA,compNA);

%compute correlation between all subjects/attachments
out = preps_multisetcca_trc(comp,'output2','single_cross','dosmooth',19);
%compute average correlation cross subjects within attachment and cross
%subject cross attachment
reps = 50;
rho = zeros(size(out.rho,2),3);
rhoshuf = zeros(size(out.rho,2),3,reps);
datn = 10;
for rep = 1:reps
    rng(rep)
    trcrep = out;
    trcrep.rho = trcrep.rho(randperm(size(trcrep.rho,1)),:);
    trc = reshape(trcrep.rho,[20,20,900]);
    dat = trc(1:10,1:10,:);
    for j = 1:size(dat,3), dat(:,:,j) = dat(:,:,j)-diag(diag(dat(:,:,j))); end
    rhoshuf(:,1,rep) = squeeze(sum(sum(dat)))./(datn.*(datn-1));
    dat = trc(11:20,1:10,:);
    for j = 1:size(dat,3), dat(:,:,j) = dat(:,:,j)-diag(diag(dat(:,:,j))); end
    rhoshuf(:,2,rep) = squeeze(sum(sum(dat)))./(datn.*(datn-1));
    rhoshuf(:,3,rep) = squeeze(mean(mean(trc(1:10,11:20,:))));
end

trc = reshape(out.rho,[20,20,900]);
dat = trc(1:10,1:10,:);
for j = 1:size(dat,3), dat(:,:,j) = dat(:,:,j)-diag(diag(dat(:,:,j))); end
rho(:,1) = squeeze(sum(sum(dat)))./(datn.*(datn-1));
dat = trc(11:20,1:10,:);
for j = 1:size(dat,3), dat(:,:,j) = dat(:,:,j)-diag(diag(dat(:,:,j))); end
rho(:,2) = squeeze(sum(sum(dat)))./(datn.*(datn-1));
rho(:,3) = squeeze(mean(mean(trc(1:10,11:20,:))));


time = out.time{1};
save(sprintf('%s/rsa/trc_parcel%d',root_dir,parcel_indx),'rho','rhoshuf','time','-v7.3')

