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
    sentences{i}  = {tmp};
    condition{i}  = stimuli_file{10}(i);
    pair_num{i}   = stimuli_file{11}(i);
    attachment{i} = stimuli_file{13}(i);
    verb_num{i}   = stimuli_file{14}(i);
    question{i}   = stimuli_file{15}(i);
    answers{i}    = stimuli_file{16}(i);
    if attachment{i} == 'NA'
    id{i}         = 100 + 10*condition{i} + 
end

stimuli = [sentences; condition; pair_num; attachment; verb_num; question; answers];

rowHeadings = {'sent_string', 'condition', 'pair_num','attachment' ...
   'verb_num', 'question','answers'};

stimuli = cell2struct(stimuli, rowHeadings,1);

%add new field for individual words
for i = 1:size(stimuli_file{1},1)%each sentence
    for j = 1:9%each word
    x{j}                        =    stimuli_file{j}(i);
    end
    rowHeadings = {'word'};
    tmp = cell2struct(x,rowHeadings,1);
    stimuli(i).words = tmp;
end
clear tmp


%save as matfile
save('/home/language/sopara/Prepositionalphrases/preps/Stimuli/preps_stimuli.m','stimuli')
