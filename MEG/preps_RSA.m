%This script loads the mscca-aligned signals from all subjects and divides
%each component into final noun-attached words and final verb-attached
%words. It then computes the between subjects correlation matrix and compares time-resolved
%correlations of noun-attached words and correlations of verb-attached
%words with correlations of between noun and verb attached words. 

root_dir    = '/project/3011210.01/MEG';
files       = dir(fullfile(root_dir,'/rsa/trc_*'));

if numel(files)<370
    for parcel_indx = 1:370
            qsubfeval('preps_execute_pipeline','preps_help_gettrc',...
                {'parcel_indx',parcel_indx},...
                'memreq',3*1024^3,'timreq',5*60);
    end
end

load atlas_subparc374_8k
pindx = 1:length(atlas.parcellationlabel);
pindx([1 2 188 189]) = []; %ignore medial wall parcels
%parcel_indx = pindx(strcmp(atlas.parcellationlabel,'R_17_B05_01'));

allrho = zeros(length(files),900,3);
allrhoshuf = zeros(length(files),900,3,50);
for parcel_indx = 1:length(files)
    load(sprintf('%s/rsa/trc_parcel%d',root_dir,parcel_indx))
    allrho(pindx(parcel_indx),:,:) = rho;
    allrhoshuf(pindx(parcel_indx),:,:,:) = rhoshuf;
end
%compute difference
diff = allrho(:,:,1) - allrho(:,:,3);
diffshuf = squeeze(allrhoshuf(:,:,1,:) - allrhoshuf(:,:,3,:));

source                = [];
source.brainordinate  = atlas;
source.label          = atlas.parcellationlabel;
source.time           = time;
source.dimord         = 'chan_time';
source.pow            = diff - mean(diffshuf,3);
cfgp                  = [];
cfgp.funparameter     = 'pow';
ft_sourcemovie(cfgp, source);
