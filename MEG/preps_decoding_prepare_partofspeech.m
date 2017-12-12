
% This script takes preprocessed data and prepares it for
% preps_decodingpipeline

load('/project/3011210.01/MEG/pilot_data_clean')

%filter

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

%extract trials Det vs Nouns
trigger_det        = [110,113,120,123,210,213,220,223]; %Determiner
trigger_noun       = [111,114,121,124,211,214,221,224]; %noun

%trigger_VA         = [218,228];
%trigger_NA         = [118,128];

[nchannels,nsamples] = size(data.trial{1});

labels             = zeros(size(data.trial,2),1);
examples           = zeros(size(data.trial,2),nchannels*nsamples);
vocab              = {'Det','Noun'};

trialsreject = [];
for trial = 1:size(data.trial,2)
    
    tmp                   = data.trial{trial};
    examples(trial,:)     = tmp(:)';
    if any(ismember(trigger_det,data.trialinfo(trial)))
        labels(trial) = 1;
    elseif any(ismember(trigger_noun,data.trialinfo(trial)))
        labels(trial) = 2;
    else
        trialsreject = [trialsreject trial];
    end
    
end
clear tmp

examples(trialsreject,:)   = [];
labels(trialsreject)    = [];

data = examples;

save(strcat('/project/3011210.01/MEG/pilot_DetvsNoun'), 'data','labels','vocab','-v7.3')   

  