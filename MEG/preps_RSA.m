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

alltrc = zeros(length(files),900,3);
for i = 1:length(files)
    load(sprintf('%s/rsa/trc_parcel%d',root_dir,i))
    alltrc(pindx(parcel_indx),:,:) = trc;
end
%compute difference
diff = alltrc(:,:,1) - alltrc(:,:,3);

    
