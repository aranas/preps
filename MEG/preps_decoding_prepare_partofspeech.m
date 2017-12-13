
% This script takes preprocessed data and prepares it for
% preps_decodingpipeline

load('/project/3011210.01/MEG/pilot_data_300')


%extract trials Det vs Nouns
trigger_det        = [110,113,120,123,210,213,220,223]; %Determiner
trigger_noun       = [111,114,121,124,211,214,221,224]; %noun
%trigger_verb       = [112,122,212,222];
%trigger_adj        = [117,127,217,227];
%trigger_prep       = [115,125,215,225];

%trigger_VA         = [218,228];
%trigger_NA         = [118,128];
%trigger_Fill       = [138];

[nchannels,nsamples] = size(data_300.trial{1});

labels             = zeros(size(data_300.trial,2),1);
examples           = zeros(size(data_300.trial,2),nchannels*nsamples);
vocab              = {'Det','Noun'};

trialsreject = [];
for trial = 1:size(data_300.trial,2)
    
    tmp                   = data_300.trial{trial};
    examples(trial,:)     = tmp(:)';
    if any(ismember(trigger_det,data_300.trialinfo(trial)))
        labels(trial) = 1;
    elseif any(ismember(trigger_noun,data_300.trialinfo(trial)))
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

  