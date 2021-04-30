function [seltrig, selpos, data, cfg] = preps_selectdata(cfg, varargin)
% This function helps select the part of the data that is relevant for the
% given analysis.
% If inputted only with one argument, the function assumes this to be
% trigger information either in form of a trigger number or a
% part-of-speech code. It will then output the corresponding trigger
% number(s) and the part of speech
% If a second argument is present, it will consider this to be the data and
% will automatically select only those parts of the data that are needed. 
% Furterh options like selecting only contentwords from the data (which is
% not coded in trigger codes for the fillers for example

% optional arguments:
% 'output_data'
% 'contentwords_only'
% 'clean_muscle'
%


data = ft_getopt(varargin, 'output_data', []);
if ~isempty(data), outdata = 1; end
content = ft_getopt(varargin, 'contentwords_only', 0);
clean = ft_getopt(varargin, 'clean_muscle', 0);

pos         = {'NA','VA','ART','NN','VVFIN','ADJA','APPR','Fill'};
trigger     = { [219,229],                     %last word Noun attached
                [119,129],                                %last word Verb attached
                [111,114,121,124,211,214,221,224,117,127,217,227], %Determiner
                [112,115,122,125,212,215,222,225,119,129,219,229],%Nouns
                [113,123,213,223],                        %Verbs
                [118,128,218,228],                        %Adjectives
                [116,126,216,226],                        %Preposition
                [31:39]};                                 %all words in filler sentences
%[40 140 240] questions

if ~iscell(cfg.seltrig)
    cfg.seltrig = {cfg.seltrig};
end
% if none are specified, take all possible triggers
if isempty(cfg.seltrig{1}) | strcmp(cfg.seltrig, 'all')
    seltrig = horzcat(trigger{:})';
    selpos = {'NA','VA','ART','NN','VVFIN','ADJA','APPR','Fill'};
    % if pos is specified, transform to trigger number
elseif isa(cfg.seltrig{1},'double')
    selpos = cell(1,length(cfg.seltrig{1}));
    for i = 1:length(trigger)
        selind = find(ismember(cfg.seltrig{1},trigger{i}));
        if ~isempty(selind), selpos(selind) = {pos{i}};end
    end
    seltrig = cfg.seltrig{1}';
    selpos = selpos';
elseif isa(cfg.seltrig{1},'char')
    selpos = {};
    for c = 1:length(cfg.seltrig)
        tmppos = cfg.seltrig{c};
        seltrig{c} = trigger{strcmp(pos,tmppos)};
        selpos = [selpos; repmat({tmppos},length(seltrig{c}),1)];
    end
    seltrig = horzcat(seltrig{:})';
else
    warning('could not collect trigger for data selection')
end
    

if outdata
    sel = ones(length(data.trialinfo),1);
    
    if clean %reject trials with muscle artifacts
        fprintf('loading artifact file for subj %s\n',cfg.subj)
        artfctfile      = fullfile('/project/3011210.01/MEG/',strcat(cfg.subj,'_muscle'));
        load(artfctfile);
        fprintf('rejecting trials with high frequency muscle noise\n')
        sel         = sel & ~ismember(data.trialinfo(:,3),noisy_trials(:,end));
        cfg.suffix      = [cfg.suffix '_cleaned'];
    end
    
    % make trial selection
    sel             = sel & ismember(data.trialinfo(:,1),seltrig);
    cfg2             = [];
    cfg2.trials      = sel;
    data         = ft_selectdata(cfg2,data);
    %fix:why does this step also slightly shift time axis?
    tn = size(data.time{1},2);
    for i = 1:length(data.trial)
        data.time{i}(1:tn) = data.time{1}(1:tn);
    end
    
        
    if content
        % not ready yet because pos in fillers are not always assigned
        % correctly 
%         load preps_stimuli
%         sent_id = unique(data.trialinfo(:,2));
%         selpos = cell(1,length(data.trialinfo));
%         for sent = 1:length(sent_id)
%             idx = ismember(data.trialinfo(:,2),sent_id(sent));
%             data.trialinfo(idx,:)
%             selpos(idx) = [stimuli(sent_id(sent)).words.pos]';
%         end        
    end
    
end

end
