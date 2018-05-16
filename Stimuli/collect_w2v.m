%This script is meant to be run only once
%It selects only those word2vec embeddings that correspond to
%words used in the experiment (saves memory & time later on)
%It then saves the selected embeddings to preps_w2v.txt which is to be used
%for all further processing

subj          = 2;

fid           = fopen(strcat('/project/3011210.01/logfiles/',num2str(subj),'_log.txt'));
format        = ['%f %f %s %f %f %f %s %f %f','%*[^\n]'];
logtxt        = textscan(fid,format,'Headerlines',6,'Delimiter','\t');

        
rowHeadings   = {'Subject','Trial','Word','Condition','PairNum'...
                'VerbNum','Attachment','Time','Duration'};
logtxt        = cell2struct(logtxt,rowHeadings,2);


% Delete questions both from words vector as well as data
ind_q           = [];

for i = 1:size(logtxt.Word,1)
    
    if logtxt.Condition(i) ~= 4
        
        ind_q   = [ind_q i];
    
    end
end


allwords        = logtxt.Word;                 %Keep all words

allwords        = allwords(ind_q);             %except questions

% clean words from sentence final period
allwords = strrep(allwords,'.','');

%extract word embeddings from pre-trained model and save to vector with
%same nrows as data.trial

fid = fopen('/project/3011210.01/word2vec/wiki.de.vec');
format = ['%s ' repmat('%f ', [1 300]),'%*[^\n]'];
wiki_embeddings_300 = textscan(fid,format);
wiki_vocab = wiki_embeddings_300(1);
wiki_vecs  = wiki_embeddings_300(2:end);
wiki_vecs  = cell2mat(wiki_vecs);
clear wiki_embeddings_300

fid = fopen('/project/3011210.01/Presentation/Stimuli/tagged_stimlist.txt');
format = '%s %s';
pos_list = textscan(fid,format,'Delimiter',' ');
pos_list{1} = strrep(pos_list{1},'.','');

unique_ind = [];
feat = [];
word = {};
pos  = {};
for i = 1:length(allwords)
        ind = find(strcmpi(wiki_vocab{1}, allwords{i}));
        if ~isempty(ind) && ~ismember(ind,unique_ind)
            unique_ind  = [unique_ind ind];
            feat        = [feat; wiki_vecs(ind,:)];
            word{end+1} = wiki_vocab{1}{ind};
            ind_pos     = find(strcmpi(pos_list{1}, allwords{i}));
            pos{end+1}  = pos_list{2}{ind_pos(1)};
        elseif isempty(ind)
            i
            allwords{i}
        end
end
w2v.feat = feat;
w2v.word = word;
w2v.pos = pos;

save('/home/language/sopara/Prepositionalphrases/preps/Stimuli/preps_w2v.mat','w2v')