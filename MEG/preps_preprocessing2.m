clear all
clc
% This script takes preprocessed data and prepares it for
% preps_decoding_pipeline
%% specify variables

subj            = 'pilot-002';

lastwordtrigger = [128,118,218,228];

%% load data and logfile for given subject

load(strcat('/project/3011210.01/MEG/',subj,'_dataclean'))
clear compds badcomp

%fid           = fopen(strcat('/project/3011210.01/logfiles/pilots/',num2str(subj),'_log.txt'));
fid                  = fopen(strcat('/project/3011210.01/logfiles/pilots/2_log.txt'));
format               = ['%f %f %s %f %f %f %s %f %f','%*[^\n]'];
logtxt               = textscan(fid,format,'Headerlines',6,'Delimiter','\t');


rowHeadings          = {'Subject','Trial','Word','Condition','PairNum'...
    'VerbNum','Attachment','Time','Duration'};
logtxt               = cell2struct(logtxt,rowHeadings,2);


%% select trials
label{1}            = 'Det';
trigger(1,:)        = [110,113,120,123,210,213,220,223,0]; %Determiner
label{2}            = 'Noun';
trigger(2,:)        = [111,114,121,124,211,214,221,224,0]; %noun
label{3}            = 'Verb';
trigger(3,:)        = [112,122,212,222,0,0,0,0,0]; %verbs
label{4}            = 'Adj';
trigger(4,:)        = [117,127,217,227,0,0,0,0,0]; %adjectives
label{5}            = 'Prep';
trigger(5,:)        = [115,125,215,225,0,0,0,0,0]; %prepositions
label{6}            = 'VA';
trigger(6,:)        = [218,228,0,0,0,0,0,0,0]; %last word of verb attached sentences
label{7}            = 'NA';
trigger(7,:)        = [118,128,0,0,0,0,0,0,0]; %last word of noun attached sentences
label{8}            = 'Filler';
trigger(8,:)        = [130,131,132,133,134,135,136,137,138]; %last word of filler

%% filter & save each selection to separate file

for i = 1:size(trigger,1)
    
    ind             = ismember(data.trialinfo(:,1),trigger(i,:));
    
    cfg             = [];
    cfg.trials      = ind;
    data_sel           = ft_selectdata(cfg,data);
    
    cfg             = [];
    cfg.continuous  = 'yes';
    cfg.lpfilter    = 'yes';
    cfg.lpfreq      = 30;
    cfg.lpfilttype  = 'firws';
    cfg.hpfilter    = 'yes';
    cfg.hpfreq      = 1;
    cfg.hpfilttype  = 'firws';
    cfg.usefftfilt  = 'yes';
    data_sel        = ft_preprocessing(cfg,data_sel);
    
    
    save(strcat('/project/3011210.01/MEG/',subj,'_',label{i}), 'data_sel','-v7.3')
    
end