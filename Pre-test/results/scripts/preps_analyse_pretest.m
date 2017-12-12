%Analyze pre-test data
%Assumes csv file with all numeric data of the following form: 
% column1 "Questioncode" = first number is either 1 (NA) or 2 (VA)
% following numbers are questioncode
% column2 "type of attachment = either 1 (noun) or 2 (verb)
% column3 "type of question" = either 1 (attachment) or 2 (plausibility)
% column4 "type of measurement" = either 1 (answer given) or 2 (time)
% column5... = one column per subject
path    = '/project/3011210.01/Limesurvey/responses/';
data    = csvread(strcat(path,'Aresponses.csv'),1); %skip header in first row
[n,m]   = size(data);
nsubj   = m-4;
%% Preparing data
%split data according to question type and measurement type
Plaus_d           = data(find(data(:,3) == 2 & data(:,4)==1),:);
Plaus_rt          = data(find(data(:,3) == 2 & data(:,4)==2),:);
Attach_d          = data(find(data(:,3) == 1 & data(:,4)==1),:);
Attach_rt         = data(find(data(:,3) == 1 & data(:,4)==2),:);
%
Nattach_d          = data(find(data(:,3) == 1 & data(:,4)==1 & data(:,2)==1),:);
Nattach_rt         = data(find(data(:,3) == 1 & data(:,4)==2 & data(:,2)==1),:);
Vattach_d          = data(find(data(:,3) == 1 & data(:,4)==1 & data(:,2)==2),:);
Vattach_rt         = data(find(data(:,3) == 1 & data(:,4)==2 & data(:,2)==2),:);
Nplaus_d           = data(find(data(:,3) == 2 & data(:,4)==1 & data(:,2)==1),:);
Vplaus_d           = data(find(data(:,3) == 2 & data(:,4)==1 & data(:,2)==2),:);
%% Exclude outliers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%unambiguous sentences
%List B sentences
unamb_code      = [189,159,174,157,160,171,162,192,198,169,181,152,154,173,188,167,180,1100,153,199,175,179,186];
%List A sentences
unamb_code      = [129,138,150,132,121,142,140,131,18,149,128,15,147,12,16,17,112,117,122,123,126,137];
ind             = find(ismember(Nattach_d(:,1),unamb_code'));
Nattach_unamb   = Nattach_d(ind,:);
%percentage of correct answers on unambiguous items per subject
sum(Nattach_unamb(:,5:end))/length(Nattach_unamb)
%% Finding biases 
%average time per subject for attachment decision
mean(Attach_rt(:,5:end))
std(Attach_rt(:,5:end))
%overall
mean(mean(Attach_rt(:,5:end)))
%average time per subject for attachment decision split between noun and
%verb attachments
figure;
errorbar([mean(Nattach_rt(:,5:end));mean(Vattach_rt(:,5:end))]',[std(Nattach_rt(:,5:end));std(Vattach_rt(:,5:end))]')
title('mean reaction time per subject')
legend('nouns','verbs')
%percentage of noun responses per participant and overall
tmp = Attach_d(:,5:end);
sum(tmp)/length(tmp)
sum(tmp(:))/length(tmp(:))
%% Assess material
%percentage correct answers per item (across subjects)
%count Noun answers divide by total answers
Nattach_correct    = sum(Nattach_d(:,5:end),2)/nsubj;
Vattach_correct    = 1-(sum(Vattach_d(:,5:end),2)/nsubj);
%items with below average accuracy
threshhold         = mean([Vattach_correct; Nattach_correct]);
items_fail         = Nattach_d(Nattach_correct < threshhold,1);
items_fail         = [items_fail; Vattach_d(Vattach_correct < threshhold,1)];
%RTs binned per hit/miss per subject
Nhits = Nattach_d(:,5:end)==1;
Vhits = Vattach_d(:,5:end)==0;
split_rts = [];
for i = 1:nsubj
    split_rts(i,1) = mean(Nattach_rt(Nhits(:,i),4+i));
    split_rts(i,2) = mean(Nattach_rt(~Nhits(:,i),4+i));
    split_rts(i,3) = mean(Vattach_rt(Vhits(:,i),4+i));
    split_rts(i,4) = mean(Vattach_rt(~Vhits(:,i),4+i));
end
figure;
plot(split_rts)
title('reaction times per subject')
legend('hits noun','misses noun','hits verb','misses verb')
%average response time per item
figure;
title('average RT per item')
plot(sort(mean(Nattach_rt(:,5:end),2)))
hold on
plot(sort(mean(Vattach_rt(:,5:end),2)))
legend('nouns','verbs')
%%Plausibility



