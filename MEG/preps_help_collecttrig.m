function [seltrig, selpos] = preps_help_collecttrig(subj, seltrig)

pos         = {'ART','NN','VVFIN','ADJA','APPR','NA','VA','Fill'};
trigger     = {[111,114,121,124,211,214,221,224,117,127,217,227], %Determiner
                [112,115,122,125,212,215,222,225,119,129,219,229],%Nouns
                [113,123,213,223],                        %Verbs
                [118,128,218,228],                        %Adjectives
                [116,126,216,226],                        %Preposition
                [219,229],                                %last word Noun attached
                [119,129],                                %last word Verb attached
                [30:39]};                                 %all words in filler sentences
                %[40 140 240] questions
if strcmp(subj,'pilot-002')
    warning('need to adjust trigger info')
end

if ~iscell(seltrig)
    seltrig = {seltrig};
end
%if none are specified, take all possible triggers
if isempty(seltrig{1})
    seltrig = {horzcat(trigger{:})};
end

if isa(seltrig{1},'double')
    selpos = cell(1,length(seltrig{1}));
    for i = 1:length(trigger)
        selind = find(ismember(seltrig{1},trigger{i}));
        if ~isempty(selind), selpos(selind) = {pos{i}};end
    end
    seltrig = seltrig{1}';
    selpos = selpos';
elseif isa(seltrig{1},'char')
    selpos = {};
    for c = 1:length(seltrig)
        tmppos = seltrig{c};
        seltrig{c} = trigger{strcmp(pos,tmppos)};
        selpos = [selpos; repmat({tmppos},length(seltrig{c}),1)];
    end
    seltrig = horzcat(seltrig{:})';
else
    warning('could not collect trigger for data selection')
end
    

