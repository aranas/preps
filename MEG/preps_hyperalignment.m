%% scritp computes and saves spatial filters found through multiset canonical correlation analysis

if ~exist('nfold',            'var'), nfold          = 5;                end
if ~exist('parcel_indx',      'var'), error('a parcel index needs to be specified');  end
if ~exist('datasuffix',       'var'), datasuffix    = '';                end

root_dir     = '/project/3011210.01/MEG/';

subjects = strsplit(sprintf('sub-%.3d ', [1:10]));
subjects = subjects(~cellfun(@isempty, subjects));

groupdata     = cell(1,numel(subjects));
subjectdata   = cell(1,numel(subjects));
subjecttiming = cell(1,numel(subjects));
%add for loop for parcels?
for s = 1:numel(subjects)
    % load in the data
    subj = subjects{s};
    channelfile     = fullfile(root_dir,sprintf('%s_dataclean%s.mat',subj, datasuffix));
    lcmvfile        = fullfile(root_dir,strcat(subj,'_preps_lcmv_parc.mat'));
    artfctfile      = fullfile(root_dir,strcat(subj,'_muscle'));
    
    load(channelfile,'data')
    load(artfctfile);
    artfct = find(ismember(data.trialinfo(:,3),noisy_trials(:,end)));
    %clean data % only select trials with words presented on screen
    cfg = [];
    cfg.artfctdef.reject = 'nan';
    cfg.artfctdef.muscle.artifact = data.sampleinfo(artfct,:);
    data = ft_rejectartifact(cfg, data);
    cfg = [];
    cfg.trials = ~ismember(data.trialinfo(:,1),[40 140 240]);
    data = ft_selectdata(cfg,data);
    
    load(lcmvfile);
    
    % convert the sensor-level data into  parcel-level data, for the
    % requested
    source_parc.filterlabel = filterlabel;
    subjectdata{s} = preps_sensor2parcel(data,source_parc,5,parcel_indx);
    %demean each trial
    for kk = 1:numel(subjectdata{s}.trial)
        tmp = subjectdata{s}.trial{kk};
        tmp = tmp - nanmean(tmp,2)*ones(1,size(tmp,2));
        subjectdata{s}.trial{kk} = tmp;
    end
    
    % sort subject-specific data in terms of sentence ID to match across
    % subjects
    [~, I] = sort(subjectdata{s}.trialinfo(:,2));
    groupdata{s} = subjectdata{s};%FIXME not really necessary to use new cell array here
    groupdata{s}.trialinfo = subjectdata{s}.trialinfo(I,:);
    groupdata{s}.sampleinfo = subjectdata{s}.sampleinfo(I,:);
    groupdata{s}.trial = subjectdata{s}.trial(I);
    groupdata{s}.time = subjectdata{s}.time(I);
    
    cfg            = [];
    cfg.method     = 'acrosschannel';
    groupdata{s} = ft_channelnormalise(cfg, groupdata{s});
end % for s of subjects

rng('default'); % reset the number generator, in order to be able to compare across parcels

tmpdata              = preps_groupdata2singlestruct(groupdata, subjects);
[W, A, rho, C, comp] = preps_multisetcca(tmpdata, 5, 4, []);

[comp, rho]          = preps_multisetcca_postprocess(comp, rho, source_parc.label{parcel_indx});

comp                 = ft_struct2single(comp);

savedir = '/project/3011210.01/MEG/mscca';

filename = fullfile(savedir, sprintf('mscca_parcel%03d%s',parcel_indx));
save(filename, 'rho', 'W', 'A', 'comp');
