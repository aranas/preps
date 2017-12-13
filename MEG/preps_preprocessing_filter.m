
load('/project/3011210.01/MEG/pilot_data_clean')

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
data           = ft_preprocessing(cfg,data_clean);


cfg                 = [];
cfg.resamplerefs    = 300;
cfg.demean          = 'yes';
cfg.baselinewindow  = [0.4 0.5];
data_300                = ft_resampledata(cfg,data);

save('/project/3011210.01/MEG/pilot_data_300', 'data_300','-v7.3')