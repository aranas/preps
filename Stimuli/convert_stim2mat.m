% This script takes the final stimulus list and converts it to a matlab
% file for easier access

fid = fopen('/project/3011210.01/Presentation/Stimuli/stimlist_final.txt');
format = [repmat('%s  ', [1 9]) '%f %f %f %s %f %s %s ','%*[^\n]'];
stimuli_file = textscan(fid,format,'Headerlines',1,'Delimiter','\t');

% combine words into sentence string
for i = 1:size(stimuli_file{1},1)%each sentence
    tmp = '';
    for j = 1:9%each word
        %convert from cell to str
        x = stimuli_file{j}(i);
        %concat strings
        if j == 1
          tmp = [tmp x{1}];
        else
        tmp = [tmp ' ' x{1}];
        end
    end
    id{i}         = i;
    sentences{i}  = {tmp};
    condition{i}  = stimuli_file{10}(i);
    pair_num{i}   = stimuli_file{11}(i);
    attachment{i} = stimuli_file{13}(i);
    verb_num{i}   = stimuli_file{14}(i);
    question{i}   = stimuli_file{15}(i);
    answers{i}    = stimuli_file{16}(i);
end

stimuli = [id; sentences; condition; pair_num; attachment; verb_num; question; answers];

rowHeadings = {'id','sent_string', 'condition', 'pair_num','attachment' ...
   'verb_num', 'question','answers'};

stimuli = cell2struct(stimuli, rowHeadings,1);

fid = fopen('/project/3011210.01/Presentation/Stimuli/tagged_stimlist.txt');
format = '%s %s';
pos_list = textscan(fid,format,'Delimiter',' ');

%add new field for individual words & get pos from tagged list
count = 1;
for i = 1:size(stimuli_file{1},1)%each sentence
    for j = 1:9%each word
    x{j}   = stimuli_file{j}(i);
    pos{j} = pos_list{2}(count);
    count = count + 1;
    end
    rowHeadings = {'word', 'pos'};
    tmp = cell2struct([x; pos],rowHeadings,1);
    stimuli(i).words = tmp;
    count = count + 1;
end
clear tmp


%%  add info from pre-test (added on May 14th 2018)
fid = fopen('/project/3011210.01/Presentation/Stimuli/pretest_values_per_item.txt');
format = ['%s %s %s %f %f','%*[^\n]']; %attachment, penultimate word(adj), final noun, accuracy & plausibility ratings across 20 pre-test subjects
pretest = textscan(fid,format,'Headerlines',1,'Delimiter','\t');

ind_VA = find(strcmp('VA',pretest{1}));
ind_NA = find(strcmp('NA',pretest{1}));

%fix some acc/plaus values are missing, need to be extracted from R
%notebook
for i = 1:length(stimuli)
    if strcmp('VA',stimuli(i).attachment);
        ind             = ind_VA(strcmp(strrep(stimuli(i).words(9).word,'.',''),pretest{3}(ind_VA)));
    else strcmp('NA',stimuli(i).attachment);
        ind             = ind_NA(strcmp(strrep(stimuli(i).words(9).word,'.',''),pretest{3}(ind_NA)));
    end
    if length(ind) > 1;
        ind = ind(strcmp(strrep(stimuli(i).words(8).word,'.',''),pretest{2}(ind)));
    end
    stimuli(i).acc  = pretest{4}(ind);
    stimuli(i).plaus = pretest{5}(ind);
end

for i = 1:length(stimuli)
stimuli(i).words(9).word = strrep(stimuli(i).words(9).word,'.','');
end


%save as matfile
save('/home/language/sopara/Prepositionalphrases/preps/Stimuli/preps_stimuli.mat','stimuli')

%% add frequency & co-occurrence info (added August 2018) based largely on wiki corpora
load('home/language/sopara/Prepositionalphrases/preps/Stimuli/preps_stimuli.mat')
%load lemmas and add to stimulus mat
fid = fopen('/project/3011210.01/Presentation/Stimuli/lemma_stimlist.txt');
format = ['%s %u %*[^\n]']; 
lemmas = textscan(fid,format,'Delimiter','.');
for i = 1:length(lemmas{1})
    ind = lemmas{2}(i);
    sent = strsplit(lemmas{1}{i});
    for w = 1:9
        stimuli(ind).words(w).lemma = {sent{w}};
    end
end
fclose(fid)
%load coocurrences
fid = fopen('/project/3011210.01/semanticP600/corrected_csvfiles/final_cooccurrences_correct.csv');
tLines = fgets(fid);
numCols = numel(strfind(tLines,','));
format = ['%s' repmat('%u', [1 numCols]),'%*[^\n]']; 
cooc = textscan(fid,format,'Delimiter',',');
vocab = cooc{1};
cooc  = cell2mat(cooc(2:end));
fclose(fid);
%load frequencies
fid = fopen('/project/3011210.01/semanticP600/wiki_sdewac_corpus_data/uniquestimuliwords_addedfrequency.txt');
format = ['%*s %u %*[^\n]']; 
freqs = textscan(fid,format,'Delimiter',' ');
freqs = freqs{1};
fclose(fid);

%for each sentence find corresponding freq and cooccurrence
for i = 1:length(stimuli)
    if stimuli(i).condition ~=3 %ignore fillers
    tmpcooc = zeros(9);
    tmpfreqs= zeros(1,9);
    for nword = 1:9
        tmp = find(strcmpi(stimuli(i).words(nword).lemma{1},vocab));
        if ~isempty(tmp)
            indword(nword) = tmp;
        end
    end
    whichwords = find(indword);
    tmpcooc(whichwords,whichwords) = cooc([indword(whichwords)],[indword(whichwords)]);
    stimuli(i).cooc = tmpcooc;
    %use freq values on diagonal instead of cooc
    tmpfreqs(whichwords) = freqs(indword(whichwords));
    tmpcooc = tmpcooc - diag(diag(tmpcooc)); 
    tmpcooc = tmpcooc + diag(tmpfreqs);
    end
end

save('/home/language/sopara/Prepositionalphrases/preps/Stimuli/preps_stimuli.mat','stimuli')







