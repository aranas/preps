
%Preprocessing (get each word with 100ms pre-onset and 500ms)
cfg                     = [];
cfg.dataset             = strcat('/project/3011210.01/raw/301121001sopara11_1200hz_20171207_01.ds');
cfg.logfile             = strcat('/project/3011210.01/logfiles/1_log.txt');
cfg.trialdef.prestim    = 0.1;
cfg.trialdef.poststim   = 0.5;
cfg.trialdef.eventtype  = 'UPPT001';
cfg.trialdef.eventvalue = [110:118,120:128,210:218,220:228]; 
new_cfg                  = ft_definetrial(cfg);


new_cfg.channel          = 'MEG';
new_cfg.continuous       = 'yes';
data                     = ft_preprocessing(new_cfg);

%   
% trig           = [111,112,114,121,122,124,211,212,214,221,222,224]; %trigger numbers corresponding to verbs & nouns only
% ind_noun       = find(ismember(data.trialinfo,[111,114,121,124,211,214,221,224]));
% ind_all        = ismember(data.trialinfo,trig);
% ind_all(ind_noun(201:end)) = 0;
%   
% cfg            = [];
% cfg.trials     = ind_all';
% datasel        = ft_selectdata(cfg,data);
  
cfg            = [];
cfg.continuous = 'yes';
cfg.lpfilter   = 'yes';
cfg.lpfreq     = 40;
cfg.lpfilttype = 'firws';
cfg.padding    = 10;
cfg.hpfilter   = 'yes';
cfg.hpfreq     = 0.5;
cfg.hpfilttype = 'firws';
cfg.usefftfilt = 'yes';
datasel           = ft_preprocessing(cfg,data);


