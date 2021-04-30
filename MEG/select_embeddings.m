%This function takes some data, selects corresponding word embeddings from
%stimulus file and outputs pruned data (only keeping those that have a word
%embedding)

function [brain_perword, model, all_pos, all_wids, all_embd] = select_embeddings(datasel,pos,doaverage)

%select word embeddings corresponding to selected trials
    load preps_stimuli
    model = [];
    all_pos = {};
    all_wids = {};
    all_embd = {};
    for i = 1:max(length(datasel)-1,1)
        wrd2trl = cell(length(datasel{i}.trialinfo),1);
        embd = [];
        PoS = [];
        wids = [];
        words = {};
        for w = 1:length(datasel{i}.trialinfo)
            wid = datasel{i}.trialinfo(w,2);
            wnum = num2str(datasel{i}.trialinfo(w,1));
            wnum = str2double(wnum(3));
            if ~isempty(stimuli(wid).words(wnum).w2v)
                if doaverage
                    if ~any(ismember(words,stimuli(wid).words(wnum).lemma))
                        words = [words stimuli(wid).words(wnum).lemma];
                        embd = [embd; stimuli(wid).words(wnum).w2v];
                        PoS = [PoS wnum];
                        wids = [wids wid];
                    end
                    id = find(ismember(words,stimuli(wid).words(wnum).lemma));
                    wrd2trl{id} = [wrd2trl{id} w];
                elseif ~doaverage
                    words = [words stimuli(wid).words(wnum).lemma];
                    embd = [embd; stimuli(wid).words(wnum).w2v];
                    PoS = [PoS wnum];
                    wids = [wids wid];
                end
            else
                embd = [embd; nan(1,size(stimuli(1).words(1).w2v,2))];
                PoS = [PoS nan];
                wids = [wids nan];
            end
        end
        all_pos{i} = PoS;
        all_wids{i} = wids;
        all_embd{i} = embd;
        %model RDM based on embeddings
        model_RDM = squareform(pdist(zscore(embd')','euclidean'));
        %figure;imagesc(model_RDM)
        model.rdm{i} = model_RDM;
        upos = uniqueStrCell(pos{i});
        model.name{i} = strcat(upos{:});
    end
    clear stimuli model_RDM
    
    %If models have different amount of words (due to missing embeddings
   [~,idxnan] = find(isnan(vertcat(all_pos{:})));
   if ~isempty(idxnan)
       for i = 1:length(model.rdm)
           model.rdm{i}(idxnan,:) = [];
           model.rdm{i}(:,idxnan) = [];
           all_embd{i}(idxnan,:) = [];
           all_pos{i}(idxnan) = [];
           all_wids{i}(idxnan) = [];
       end    
   end
    
    nword = length(model.rdm{1});
    [nparc, nt] = size(datasel{end}.trial{1});
    if doaverage
        % Average over repeated words
        brain_perword = zeros(nword,nparc,nt);
        for w = 1:length(words)
            trlsel = datasel{end}.trial(wrd2trl{w});
            brain_perword(w,:,:) = mean(cat(3,trlsel{:}),3);
        end
        clear trlsel
    else
        brain_perword = permute(cat(3,datasel{end}.trial{:}),[3 1 2]);
        brain_perword(idxnan,:,:) = [];
    end
    
%     if 1%model_attach
%         attach = datasel{end}.trialinfo(:,1);
%         attach(idxnan) = [];
%         embd = all_embd{1};
%         embd(ismember(attach,[219,229]),:) = all_embd{2}(ismember(attach,[219,229]),:);
%         
%         model.rdm{end+1} = squareform(pdist(zscore(embd')','euclidean'));
%         model.name{end+1} = 'attachment_relabeled';
%     end