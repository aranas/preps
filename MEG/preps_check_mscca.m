

% this script is for testing the results of the mscca
if ~exist('do_premscca',        'var'), do_premscca = false;            end
if ~exist('parcel_indx',        'var'), parcel_indx = [];            end

root_dir    = '/project/3011210.01/MEG';
files       = dir(fullfile(root_dir,'mscca'));
files       = files(3:end);

load atlas_subparc374_8k

pindx = 1:length(atlas.parcellationlabel);
pindx([1 2 188 189]) = []; %ignore medial wall parcels
%parcel_indx = pindx(strcmp(atlas.parcellationlabel,'R_17_B05_01'));

% alltrc = zeros(length(files),301);
% alltrcVA = zeros(length(files),900);
% alltrcNA = zeros(length(files),900);
%for parcel_indx = 1:length(files)
    
    load(fullfile(root_dir,'mscca',files(parcel_indx).name))
    
    %trc = preps_multisetcca_trc(comp,'output','Z_scaled');
    
    %figure;plot(trc.time{1},trc.rho)
    %hold on;
    
    %% compare with lcmv principal component
    if do_premscca
        subjects = strsplit(sprintf('sub-%.3d ', [1:10]));
        subjects = subjects(~cellfun(@isempty, subjects));
        datasuffix = '';
        for s = 1:numel(subjects)
            s
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
            
            load(lcmvfile,'source_parc','filterlabel');
            
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
        %select 1st principal component
        tmpdata     = preps_groupdata2singlestruct(groupdata, subjects);
        for i = 1:length(tmpdata.trial)
            tmpdata.trial{i} = tmpdata.trial{i}([1 6 11 16 21 26 31 36 41 46],:);
        end
        tmpdata.label = tmpdata.label([1 6 11 16 21 26 31 36 41 46]);

        trcpre      = preps_multisetcca_trc(tmpdata,'output','Z_scaled');
        hold on
        plot(trcpre.time{1},trcpre.rho)
    end
    %
    % legend('pre-mscca 1st pc','on first mscca component')
    % open preps_plotting
    % set(gca,'FontSize',25)
    % set(legend,'Location','best')
    % xlabel(sprintf('time in ms '))
    % ylabel('average between-subject correlation')
    % title('Time-resolved correlation across all 10 subjects for visual parcel')
    % save_file = '/project/3011210.01/MEG/figures/trc_R17_prepostmscca'
    % export_fig(save_file,'-png');
    
    %separate final words according to VA/NA
    %repeat preps_multisetcca_trc for each condition individually
    %permute trials across conditions
    cfg = [];
    cfg.trials = ismember(comp.trialinfo(:,1),[129,119]);
    compVA = ft_selectdata(cfg,comp);
    cfg.trials = ismember(comp.trialinfo(:,1),[229,219]);
    compNA = ft_selectdata(cfg,comp);
    
    trcVA = preps_multisetcca_trc(compVA,'output','Z_scaled','dosmooth',19);
    trcNA = preps_multisetcca_trc(compNA,'output','Z_scaled','dosmooth',19);
    
    %alltrc(pindx(parcel_indx),:) = trc.rho;
%     alltrcVA(pindx(parcel_indx),:) = trcVA.rho;
%     alltrcNA(pindx(parcel_indx),:) = trcNA.rho;

% rng('default')
% tmptrc = zeros(500,10,10,900);
% idVA = unique(compVA.trialinfo(:,2));
% idNA = unique(compNA.trialinfo(:,2));
% %fix some irregularities in IDs
% idNA = [idNA(1:7);175;idNA(8:end)];
% idNA = [idNA(1:25);176;idNA(26:end)];
% idNA(89:90) = [];
%
% comprand = compVA;
% for p = 1:500
%     p
%     for n = 1:length(comprand.trial)
%         indx(p,n,:) = randi(2,1,10);
%         idxVA = find(indx(p,n,:)==1);
%         tmpidx = idNA(find(idVA==comprand.trialinfo(n,2)));
%         comprand.trial{n}(find(indx(p,n,:)==2),:) = compNA.trial{find(compNA.trialinfo(:,2)==tmpidx)}(find(indx(p,n,:)==2),:);
%     end
%     tmptrc(p,:,:,:) = preps_multisetcca_trc(comprand);
% end

