clear all
clc

%Preprocessing (get each word with 100ms pre-onset and 500ms)
cfg                     = [];
cfg.dataset             = strcat('/project/3011210.01/raw/301121001sopara11_1200hz_20171207_01.ds');
cfg.logfile             = strcat('/project/3011210.01/logfiles/1_log.txt');
cfg.trialdef.prestim    = 0.5;
cfg.trialdef.poststim   = 1;
cfg.trialdef.eventtype  = 'UPPT001';
cfg.trialdef.eventvalue = [110:118,120:128,210:218,220:228,130:138]; 
new_cfg                  = ft_definetrial(cfg);


new_cfg.channel          = {'MEG', 'EEG'};
new_cfg.continuous       = 'yes';
data                     = ft_preprocessing(new_cfg);


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
    cfg = [];
    cfg.component = [1:30];       % specify the component(s) that should be plotted
    cfg.layout    = 'CTF275.lay'; % specify the layout file that should be used for plotting
    cfg.comment   = 'no';
    ft_topoplotIC(cfg, compds)

    fh = figure();
    cfg = [];
    %cfg.channel = [1:30]; % components to be plotted
    cfg.channel     = unique(j((end-2):end, :));
	cfg.compscale   = 'local';
    cfg.viewmode = 'component';
    cfg.layout = 'CTF275.lay'; % specify the layout file that should be used for plotting
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
    cfg = [];
    cfg.channel = 'MEG';
    cfg.component = badcomp; % to be removed component(s)
    data_clean = ft_rejectcomponent(cfg, comp, data);
    
    % Remove non-MEG channels, because ft_rejectcomponent adds those again
    cfg         = [];
    cfg.channel = 'MEG';
    data_clean   = ft_preprocessing(cfg, data_clean);
  
    save('/project/3011210.01/MEG/pilot_data_clean', 'data_clean','compds','badcompds','-v7.3')
      