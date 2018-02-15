
% This script takes preprocessed data and prepares it for
% preps_decodingpipeline

load('/project/3011210.01/MEG/p2_data_clean')

%% select trials
trigger_det     = [110,113,120,123,210,213,220,223]; %Determiner
trigger_noun    = [111,114,121,124,211,214,221,224]; %noun
trigger_verb    = [112,122,212,222]; %verbs
trigger_adj     = [117,127,217,227]; %adjectives
trigger_prep    = [115,125,215,225]; %prepositions

trigger_VA      = [218,228]; %last word of verb attached sentences
trigger_NA      = [118,128]; %last word of noun attached sentences
trigger_Fill    = [138]; %last word of filler

ind             = ismember(data.trialinfo,[trigger_noun,trigger_adj]);

cfg             = [];
cfg.trials      = ind;
data            = ft_selectdata(cfg,data);

%%  filter

cfg            = [];
cfg.continuous = 'yes';
cfg.lpfilter   = 'yes';
cfg.lpfreq     = 30;
cfg.lpfilttype = 'firws';
%cfg.padding    = 10; 
% cfg.hpfilter   = 'yes';
% cfg.hpfreq     = 1;
% cfg.hpfilttype = 'firws';
cfg.usefftfilt = 'yes';
data           = ft_preprocessing(cfg,data);

%%
[nchannels,nsamples] = size(data.trial{1});

labels              = zeros(size(data.trial,2),1);
vocab               = {'Noun','Adjective'};

for trial = 1:size(data.trial,2)
    if any(ismember(trigger_noun,data.trialinfo(trial)))
        labels(trial)   = 1;
    elseif any(ismember(trigger_adj,data.trialinfo(trial)))
        labels(trial)   = 2;
    end
end


save(strcat('/project/3011210.01/MEG/p2_NounvsAdj_nolof'), 'data','labels','vocab','-v7.3')

