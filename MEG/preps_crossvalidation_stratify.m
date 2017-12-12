% This function is not working!!!
%This function creates a cross-validation scheme 
%It gets as input how many examples should be left out (nleaveout), how
%many observations there are (numobs) as well as the column vector with
%labels per observations.
%It ouputs two numobsx1 column vectors with numbers 1:nfolds signaling
%which observations will be left out (testind) per
%fold

function [testind,nfolds] = preps_crossvalidation_stratify(nleaveout,labels,numobs) 

% Find frequency of each label in observations
[classLabels,~,counts]  = unique(labels);
nclass                  = hist(counts,length(classLabels));
classPriors             = nclass / sum(nclass);

% Find all possible numbers to leave out
range =1:numobs;

temp         = range.*classPriors(2);
ntest               = temp(find(floor(temp) == temp))



fprintf('creating cross-validation scheme');
fprintf('\n');

ntrain          = nobs-ntest;
nperlabeltest   = ntest*classPriors;  % these must be integers, of course
nperlabeltrain  = ntrain*classPriors;
nfolds          = floor(nobs/sum(nperlabeltest));
    
fprintf('number of examples to leave out for testing: ');
fprintf('\n');
for i = 1:size(nperlabeltest,2)
    fprintf('for Label %d: %d',i,nperlabeltest(i));
    fprintf('\n');
end
fprintf('number of folds: %d',nfolds);
fprintf('\n');


testind   = zeros(numobs,1);

for i = 1:length(classLabels) % create permedLabIdxs{i} as random perm of examples with ith label
    ind_label          = find(labels==classLabels(i));
	tmp_ind            = crossvalind('Kfold',length(ind_label),floor(length(ind_label)/nperlabeltest(1)));
    testind(ind_label) = tmp_ind;
end

        
