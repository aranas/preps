% This script extracts information from logfile about which sentence/word a
% trial belongs to. Where each trial corresponds to one trigger, thus one
% word.
% cleaned data contains all trials except questions
subj          = 2;

load(strcat('/project/3011210.01/MEG/p',num2str(subj),'_data_clean_alltrials'))
clear compds badcomp

fid           = fopen(strcat('/project/3011210.01/logfiles/',num2str(subj),'_log.txt'));
format        = ['%f %f %s %f %f %f %s %f %f','%*[^\n]'];
logtxt        = textscan(fid,format,'Headerlines',6,'Delimiter','\t');

        
rowHeadings   = {'Subject','Trial','Word','Condition','PairNum'...
                'VerbNum','Attachment','Time','Duration'};
logtxt        = cell2struct(logtxt,rowHeadings,2);


% See if words from logfile correspond to amount of triggers:
if size(logtxt.Word,1) ~= size(data.trialinfo,1)
    warning('logfile does not correspond to trigger numbers!!!')
end

% Delete questions both from words vector as well as data
ind_q           = [];

for i = 1:size(logtxt.Word,1)
    
    if logtxt.Condition(i) ~= 4
        
        ind_q   = [ind_q i];
    
    end
end


allwords        = logtxt.Word;                 %Keep all words

allwords        = allwords(ind_q);             %except questions

data            = ft_selectdata(data,'rpt',ind_q);%keep all trial that are not questions

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

feat = zeros(length(allwords),size(wiki_vecs,2));
ind_match = [];
for i = 1:length(allwords)
        ind = find(strcmpi(wiki_vocab{1}, allwords{i}));
        if ~isempty(ind)
            feat(i,:) = wiki_vecs(ind,:);
            ind_match = [ind_match i];
        end
end
% many nonmatches are due to spelling errors (should be corrected in
% stimulus set/logfiles
allwords = allwords(ind_match); 
feat     = feat(ind_match,:);
data     = ft_selectdata(data,'rpt',ind_match);%keep all trial that are not questions


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


save(strcat('/project/3011210.01/MEG/p2_alltrials_w2v'), 'data','allwords','feat','-v7.3')









%Find indices that belong to questions
%filter out questions with wrong trigger numbering
%question triggers are 119,129,219,229 & some of the 130:
%indx_130 = find(new_cfg.trl(:,4)==130);
%indx_130 = indx_130(new_cfg.trl(indx_130+1,4)~=131);
%indx_130 = indx_130(new_cfg.trl(indx_130-1,4)==138);