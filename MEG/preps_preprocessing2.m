clear all
clc
% This script takes preprocessed data and prepares it for
% preps_decoding_pipeline
subj = 'p2';
lastwordtrigger = [128,118,218,228];
%% load data and logfile for given subject

load(strcat('/project/3011210.01/MEG/',subj,'_data_clean_alltrials'))
clear compds badcomp

%fid           = fopen(strcat('/project/3011210.01/logfiles/pilots/',num2str(subj),'_log.txt'));
fid                  = fopen(strcat('/project/3011210.01/logfiles/pilots/2_log.txt'));
format               = ['%f %f %s %f %f %f %s %f %f','%*[^\n]'];
logtxt               = textscan(fid,format,'Headerlines',6,'Delimiter','\t');


rowHeadings          = {'Subject','Trial','Word','Condition','PairNum'...
    'VerbNum','Attachment','Time','Duration'};
logtxt               = cell2struct(logtxt,rowHeadings,2);

%% Take care of unwanted triggers
if strcmp(subj,'p2')
    ind              = ~ismember(data.trialinfo,192);
    cfg              = [];
    cfg.trials       = ind;
    data             = ft_selectdata(cfg,data);
end

%% assign FastText word embeddings vector


% See if words from logfile correspond to amount of triggers:
if size(logtxt.Word,1) ~= size(data.trialinfo,1)
    warning('logfile does not correspond to trigger numbers!!!')
end

% Delete questions both from words vector as well as data
ind_q                = [];

for i = 1:size(logtxt.Word,1)
    if logtxt.Condition(i) ~= 4
        ind_q        = [ind_q i];
    end
end


allwords             = logtxt.Word;                 %Keep all words
allwords             = allwords(ind_q);             %except questions
data                 = ft_selectdata(data,'rpt',ind_q);%keep all trial that are not questions

% clean words from sentence final period
allwords             = strrep(allwords,'.','');
clear logtxt ind_q
%extract word embeddings from pre-trained model and save to vector with
%same nrows as data.trial

load('/home/language/sopara/Prepositionalphrases/preps/Stimuli/preps_w2v.mat')

allfeat              = zeros(length(allwords),size(w2v.feat,2));
allpos               = cell(length(allwords),1);
ind_match            = [];
for i = 1:length(allwords)
    ind              = find(strcmpi(w2v.word, allwords{i}));
    if ~isempty(ind)
        allfeat(i,:) = w2v.feat(ind,:);
        allpos{i}    = w2v.pos{ind};
    else
        allwords{i}
    end
end
clear w2v
%% assign pre-test information
load('/home/language/sopara/Prepositionalphrases/preps/Stimuli/preps_stimuli.mat')

indices = find(ismember(data.trialinfo,lastwordtrigger)); %loop through sentences based on last word match with stimulus file

for i = 1:length(stimuli)
    ind2            = [];
    ind             = find(strcmp(stimuli(i).words(9).word,allwords));
    if strcmp(stimuli(i).attachment,'NA');
        [~,ind2]    = min(data.trialinfo(ind));
    elseif strcmp(stimuli(i).attachment,'VA');
        [~,ind2]    = max(data.trialinfo(ind));
    end
    if ~isempty(ind2) && ~isempty(stimuli(i).acc)
        for j = 0:8
            data.trialinfo(ind(ind2)-j,2) = stimuli(i).acc;
            data.trialinfo(ind(ind2)-j,3) = stimuli(i).plaus;
        end
    end
end
clear stimuli

%% select trials
label{1}            = 'Det';
trigger(1,:)        = [111,114,121,124,211,214,221,224,0]; %Determiner
label{2}            = 'Noun';
trigger(2,:)        = [112,115,122,125,212,215,222,225,0]; %noun
label{3}            = 'Verb';
trigger(3,:)        = [113,123,213,223,0,0,0,0,0]; %verbs
label{4}            = 'Adj';
trigger(4,:)        = [118,128,218,228,0,0,0,0,0]; %adjectives
label{5}            = 'Prep';
trigger(5,:)        = [116,126,216,226,0,0,0,0,0]; %prepositions
label{6}            = 'VA';
trigger(6,:)        = [218,228,0,0,0,0,0,0,0]; %last word of verb attached sentences
label{7}            = 'NA';
trigger(7,:)        = [118,128,0,0,0,0,0,0,0]; %last word of noun attached sentences
label{8}            = 'Filler';
trigger(8,:)        = [130,131,132,133,134,135,136,137,138]; %last word of filler

%% filter & save each selection to separate file

for i = 7:size(trigger,2)
    
    ind             = ismember(data.trialinfo(:,1),trigger(i,:));
    
    cfg             = [];
    cfg.trials      = ind;
    data_sel           = ft_selectdata(cfg,data);
    
    words           = allwords(ind);
    pos             = allpos(ind);
    feat            = allfeat(ind,:);

    cfg             = [];
    cfg.continuous  = 'yes';
    cfg.lpfilter    = 'yes';
    cfg.lpfreq      = 30;
    cfg.lpfilttype  = 'firws';
    cfg.hpfilter    = 'yes';
    cfg.hpfreq      = 0.5;
    cfg.hpfilttype  = 'firws';
    cfg.usefftfilt  = 'yes';
    data_sel        = ft_preprocessing(cfg,data_sel);
    
    
    save(strcat('/project/3011210.01/MEG/',subj,'_',label{i}), 'data_sel','words','pos','feat','-v7.3')
    
end