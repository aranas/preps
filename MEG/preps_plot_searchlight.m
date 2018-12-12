function preps_plot_searchlight(subj,filename)

root_dir    = '/project/3011210.01/MEG/Classification';
files       = dir(fullfile(root_dir,subj,filename));

%check if all parcels have been computed
allp  = {files(:).name};
for parcel_indx = 1:370
    if isempty(cell2mat(strfind(allp,sprintf('parcel%03d',parcel_indx))))
        qsubfeval('preps_execute_pipeline','preps_decoding',{'subj',subj},{'classifier','preps_naivebayes'},{'mode','normal'},{'dattype','lcmv'},{'seltrig',{'NA', 'VA'}},{'clean_muscle',1},{'dopca',0 },{'toverlap',0.8},{'datasuffix','_lp01'},{'numfeat','all'},{'parcel_indx',parcel_indx},{'resample',1},'memreq',8*1024^3,'timreq',5*60*60,'batchid',sprintf('preps_searchlight_NAVA_%s_parcel%03d',subj,parcel_indx));
    end
end

if strcmp(subj,'group') %combine results to get group average
    subjects = strsplit(sprintf('sub-%.3d ', [1:10]));
    subjects = subjects(~cellfun(@isempty, subjects));
    
    for s = 1:10
        subj = subjects{s};
        file = dir(fullfile(root_dir,subj,'lcmv/searchlight/*NAVA_allparcels.mat'));
        file
        load(fullfile(file.folder,file.name))
        allacc(s,:,:,:) = acc;
        allaccshuf(s,:,:,:) = accshuf;
    end
    taxis = round(cfgcv.time(1:nearest(cfgcv.time,cfgcv.time(end)-cfgcv.twidth/2)));
    clear acc accshuf
    allacc(:,[1 2 188 189],:,:) = nan;
    allaccshuf(:,[1 2 188 189],:,:) = nan;
    load atlas_subparc374_8k
    load(fullfile('/project/3011210.01/anatomy',subj,strcat(subj,'_sourcemodel.mat')));
    T                     = squeeze(nanmean(nanmean(allacc,4),1));
    source                = [];
    source.brainordinate  = atlas;
    source.brainordinate.pos = sourcemodel.pos;
    source.label          = atlas.parcellationlabel;
    source.time           = taxis;
    source.dimord         = 'chan_time';
    source.pow            = T;
    source.mask           = double(T>0.5);
    
    source.pow            = source.pow.*source.mask;
    cfgp                  = [];
    cfgp.funparameter     = 'pow';
    %cfgp.maskparameter    = 'mask';
    ft_sourcemovie(cfgp, source);
    
else
    
    load atlas_subparc374_8k
    pindx = 1:length(atlas.parcellationlabel);
    pindx([1 2 188 189]) = []; %ignore medial wall parcels
    
    for parcel_indx = 1:length(files)
        parcel_indx
        load(fullfile(files(parcel_indx).folder,files(parcel_indx).name));
        
        stat        = cell2mat(stat);
        stat        = struct2cell(stat);
        statshuf    = cell2mat(statshuf);
        statshuf    = struct2cell(statshuf);
        indxparc    = strfind(files(parcel_indx).name,'parcel');
        indx        = pindx(str2double(files(parcel_indx).name(indxparc+6:indxparc+8)));
        
        acc(indx,:,:)       = cell2mat(squeeze(stat(1,:,:)));
        accshuf(indx,:,:)   = cell2mat(squeeze(statshuf(1,:,:)));
        
        if parcel_indx==1,
            acc(374,end)=0;
            accshuf(374,end,end)=0;
        end
        
    end
    taxis = round(cfgcv.time(1:nearest(cfgcv.time,cfgcv.time(end)-cfgcv.twidth/2))*1000);
    cfgcv.time = taxis;
    indsuffix = strfind(files(1).name,'parcel')-1;
    filename = fullfile(files(1).folder,[files(1).name(1:indsuffix),'allparcels']);
    save(filename,'acc','accshuf','cfgcv');
    acc(acc==0) = nan;
    accshuf(accshuf==0) = nan;
    
    load(fullfile('/project/3011210.01/anatomy',subj,strcat(subj,'_sourcemodel.mat')));
    source                = [];
    source.brainordinate  = atlas;
    source.brainordinate.pos = sourcemodel.pos;
    source.label          = atlas.parcellationlabel;
    source.time           = taxis;
    source.dimord         = 'chan_time';
    source.pow            = nanmean(acc,3);
    source.mask           = double(mean(acc,3)>0.5);
    
    source.pow            = source.pow.*source.mask;
    
    %plot both hemispheres simultaneously
    pos = sourcemodel.pos;
    n = 7842;
    pos(1:n,1) = pos(1:n,1)+150;
    pos(1:n,1:2) = -pos(1:n,1:2);
    source.brainordinate.pos = pos;
    
    cfgp                  = [];
    cfgp.funparameter     = 'pow';
    cfgp.time             = taxis;
    %cfgp.maskparameter    = 'mask';
    ft_sourcemovie(cfgp, source);
    
end
end