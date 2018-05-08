clear all
clc
% This script takes preprocessed data and prepares it for
% preps_decodingpipeline

load('/project/3011210.01/MEG/p5_data_clean_alltrials')

%% select trials
trigger_det     = [111,114,121,124,211,214,221,224]; %Determiner
trigger_noun    = [112,115,122,125,212,215,222,225]; %noun
trigger_verb    = [113,123,213,223]; %verbs
trigger_adj     = [118,128,218,228]; %adjectives
trigger_prep    = [116,126,216,226]; %prepositions

trigger_VA      = [219,229]; %last word of verb attached sentences
trigger_NA      = [119,129]; %last word of noun attached sentences
trigger_Fill    = [139]; %last word of filler

ind             = ismember(data.trialinfo,[trigger_verb,115,125,215,225,trigger_VA,trigger_NA]);

cfg             = [];
cfg.trials      = ind;
data           = ft_selectdata(cfg,data);

%%  filter

cfg            = [];
cfg.continuous = 'yes';
cfg.lpfilter   = 'yes';
cfg.lpfreq     = 30;
cfg.lpfilttype = 'firws';
cfg.padding    = 10; 
cfg.hpfilter   = 'yes';
cfg.hpfreq     = 0.5;
cfg.hpfilttype = 'firws';
cfg.usefftfilt = 'yes';
data           = ft_preprocessing(cfg,data);

%%
[nchannels,nsamples] = size(data.trial{1});

labels              = zeros(size(data.trial,2),1);
vocab               = {'Verb','Noun'};

for trial = 1:size(data.trial,2)
    if any(ismember(trigger_verb,data.trialinfo(trial)))
        labels(trial)   = 1;
    elseif any(ismember([115,125,215,225],data.trialinfo(trial)))
        labels(trial)   = 2;
    elseif any(ismember(trigger_VA,data.trialinfo(trial)))
        labels(trial)   = 1;
    elseif any(ismember(trigger_NA,data.trialinfo(trial)))
        labels(trial)   = 2;
    end
end


save(strcat('/project/3011210.01/MEG/p5_VerbvsNoun_generalize_05hp'), 'data','labels','vocab','-v7.3')

