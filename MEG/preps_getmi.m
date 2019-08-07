function [midiff,data] = preps_getmi(data)

totalcount = 2343449936;                                    %total amount of words in dictionary
load('/project/3011210.01/semanticP600/preps_stimuli.mat')  %load stimulus matrix to identify words
stimuli = stimuli([stimuli.condition]~=3);
%threshold for how outlier frequency counts as two standard deviations away
%from mean.
freqs = [];
for i = 1:length(stimuli)
    id = stimuli(i).id;
    freqs = [freqs;diag(stimuli(id).cooc)];
end
freqs(freqs==0) = [];
thresh = mean(log(freqs))-(2*std(log(freqs)));
%% Extract difference in coocurrence values as independent variable
%For each trial
% norm_cooc1 = zeros(1,length(data.trial));
% norm_cooc2 = zeros(1,length(data.trial));
mi_cooc1 = zeros(1,length(data.trialinfo));
mi_cooc2 = zeros(1,length(data.trialinfo));
for t = 1:length(data.trialinfo)
    id = data.trialinfo(t,2);
    %1. 3rd word (verb) & 9th word(noun) 2. 5th word (noun) & 9th word (noun)
    %extract cooc and compute mutual information
    mi_cooc1(t)  = log2((stimuli(id).cooc(3,9)*totalcount)/(stimuli(id).cooc(3,3)*stimuli(id).cooc(9,9)));
    mi_cooc2(t)  = log2((stimuli(id).cooc(5,9)*totalcount)/(stimuli(id).cooc(5,5)*stimuli(id).cooc(9,9)));
    % set MI to zero, if cooccurrence is zero and individual freq are
    % above threshold
    if isinf(mi_cooc1(t)) || isinf(mi_cooc2(t))
        if log(stimuli(id).cooc(9,9)) > thresh
            if log(stimuli(id).cooc(3,3)) > thresh
                mi_cooc1(t) = 0;
            end
            if log(stimuli(id).cooc(5,5)) > thresh
                mi_cooc2(t) = 0;
            end
        end
    end
end
%compute difference in coocurrence
midiff   = abs(mi_cooc1(:) - mi_cooc2(:));

% Get rid of trials with NaNs or Infs in independent variable
non_num     = isnan(midiff) + isinf(midiff);
midiff    = midiff(~non_num);

cfg         = [];
cfg.trials  = ~non_num;
data     = ft_selectdata(cfg,data);
