clear all
clc

%% Specify variables
if ~exist('subj',        'var'), subj         = 'sub-009';  end
if ~exist('trigger',     'var'), trigger      = [140,240,40,31:39,110:119,120:129,210:219,220:229];  end
if strcmp(subj,'pilot-002')
    trigger = [130:139,110:119,120:129,210:219,220:229];
end
if ~exist('trigger_last','var'), trigger_last = [119,129,219,229];  end
if ~exist('root_dir',    'var'), root_dir     = '/project/3011210.01/';  end

%% epoch data from raw files
file                    = dir(strcat('/project/3011210.01/raw/',subj,'/ses-meg01/meg/'));

cfg                     = [];
cfg.dataset             = strcat(root_dir,'raw/',subj,'/ses-meg01/meg/',file(3).name);
cfg.logfile             = strcat(root_dir,'raw/',subj,'/ses-meg01/beh/',subj, '_log.txt');
cfg.trialdef.prestim    = 0.2;
cfg.trialdef.poststim   = 2.8;
cfg.trialdef.eventtype  = 'UPPT001';
cfg.trialdef.eventvalue = trigger;
cfg.trialfun            = 'ft_trialfun_preps';
new_cfg                 = ft_definetrial(cfg);

%% before filtering data, identify trials with muscle activity
new_cfg.continuous      = 'yes';
new_cfg.demean          = 'yes';
%new_cfg.dftfilter       = 'yes';    %get rid of line noise
data                    = ft_preprocessing(new_cfg);
data.trialinfo(:,3)     = 1:length(data.trialinfo); %add numbering to trials
% redefine to keep longer time windows only for sentence-final words
ind_lastword            = find(ismember(data.trialinfo(:,1),trigger_last));
toi                     = repmat([data.time{1}(1) 0.8],length(data.trial),1);
toi(ind_lastword,:)     = repmat([data.time{1}(1) data.time{1}(end)],length(ind_lastword),1);

cfg                     = [];
cfg.channel             = 'MEG';
cfg.latency             = toi; %adapt rejectvisual_summary
cfg.preproc.hpfilter    = 'yes';
cfg.preproc.hpfreq      = 100;
cfg.metric              = 'zvalue';
tmp_data_muscle         = ft_rejectvisual(cfg, tmpdata);
artfct                  = setdiff(data.trialinfo(:,3),tmp_data_muscle.trialinfo(:,3));
noisy_trials            = data.trialinfo(artfct,:);
% cfg = [];
% cfg.channel = 'MEG';
% ft_databrowser(cfg,tmpdata)

save(strcat(root_dir,'MEG/',subj,'_muscle'),'noisy_trials','-v7.3')
clear tmpdata tmp_data_muscle noisy_trials
%% filter
new_cfg.channel             = {'MEG', 'EEG'};
new_cfg.lpfilter            = 'yes';
new_cfg.lpfreq              = 40;
new_cfg.lpfilttype          = 'firws';
new_cfg.hpfilter            = 'yes';
new_cfg.hpfreq              = 0.1;
new_cfg.hpfilttype          = 'firws';
new_cfg.usefftfilt          = 'yes';
new_cfg.padding             = 10;
data                        = ft_preprocessing(new_cfg);

data.trialinfo(:,3)     = 1:length(data.trialinfo); %add numbering to trials

cfg                     = [];
cfg.toilim              = toi;
cfg.trials              = 'all';
data                    = ft_redefinetrial(cfg,data);
% 
% cfg = [];
% cfg.channel = 'MEG';
% ft_databrowser(cfg,data)

%% compute ICA on data to remove ECG/EOG artifacts
%downsampling data to 300 Hz for ICA analysis

cfg                 = [];
cfg.channel         = 'MEG';
cfg.resamplefs      = 300;
cfg.detrend         = 'no';
data_resamp         = ft_resampledata(cfg, data);


% ICA analysis (limit to 100 iterations)
cfg                 = [];
cfg.channel         = 'MEG';
cfg.method          = 'runica';
cfg.runica.maxsteps = 100;
compds              = ft_componentanalysis(cfg, data_resamp);

% perform same ICA on original data (not downsampled)
cfg                 = [];
cfgcomp.channel     = 'MEG';
cfg.unmixing        = compds.unmixing;
cfg.topolabel       = compds.topolabel;
comp                = ft_componentanalysis(cfg,data);


%% Correlate ICs to EOGh, EOGv and ECG and store which components to remove

EEG = cell2mat(data.trial);
EEG = EEG(ismember(data.label, {'EEG057', 'EEG058','EEG059'}), :); %or 'EEG059'

% Correlate ICs to EEG and print the ICs that have highest correlations
r = corr(cell2mat(comp.trial)', EEG');
[~, j] = sort(abs(r));
[r(j((end-3):end, 1), 1), r(j((end-3):end, 2), 2), r(j((end-3):end, 3), 3), j((end-3):end, :)]


%% STEP 3.2
% inspect ICA component and marks those related to ECG/EOG artifacts as bad
cfg           = [];
cfg.component = [1:30];       % specify the component(s) that should be plotted
cfg.layout    = 'CTF275.lay'; % specify the layout file that should be used for plotting
cfg.comment   = 'no';
ft_topoplotIC(cfg, compds)

fh            = figure();
cfg           = [];
cfg.channel   = unique(j((end-2):end, :));
cfg.compscale = 'local';
cfg.viewmode  = 'component';
cfg.layout    = 'CTF275.lay'; % specify the layout file that should be used for plotting
ft_databrowser(cfg, compds)

waitfor(fh) % wait until we have closed the window for entring bad components
%% STEP 3.2
% save bad ICA components

badcomp='tmp';
promptUser=true;
while promptUser
    prompt=inputdlg('List of ICA components to remove','Output File',1,{'tmp'});
    if isempty(prompt)
        disp('Cancel experiment ...');
        return;
    else
        badcomp= str2num(prompt{1});
    end
    
    if badcomp
        promptUser = false;
    end
end
%% STEP 3.3
% remove components related to ECG and EOG artifacts and backproject the data
cfg                 = [];
cfg.channel         = 'MEG';
cfg.component       = badcomp; % to be removed component(s)
data                = ft_rejectcomponent(cfg, comp, data);

cfg                 = [];
cfg.resamplefs      = 300;
cfg.detrend         = 'no';  % not good for evoked data
cfg.demean          = 'no';
cfg.trials          = 'all';
data                = ft_resampledata(cfg, data);

% Remove non-MEG channels, because ft_rejectcomponent adds those again
cfg                 = [];
cfg.channel         = 'MEG';
cfg.detrend         = 'no';
cfg.demean          = 'yes';
cfg.baselinewindow  = [-inf 0];
data                = ft_preprocessing(cfg, data);


save(strcat(root_dir,'MEG/',subj,'_dataclean'),'data','compds','badcomp','-v7.3')
%missing = [ 16 65 163 199 212   234   266   338   361   366];