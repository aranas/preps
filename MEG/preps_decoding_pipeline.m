% This script trains a classifier/model to some data with the
% option to apply a cross-validation scheme.
% This script assumes two variables:
% data      = mxn matrix, which is a feature matrix with m trials/observations 
%               and n features.
% labels    = mx1 column vector containing the labels for each trial/observation.
% 
addpath('/home/language/sopara/Prepositionalphrases/Scripts_tom/GNB')
addpath('/home/language/sopara/Prepositionalphrases/preps/MEG/')
%% Script variables (to be specified before running script)
if ~exist('pathname', 'var')
    error('a pathname needs to be provided');
end

if ~exist('rootdir', 'var')
    rootdir = '/project/3011210.01/';
end

if ~exist('dofolds','var'),         dofolds         = 1; end
if ~exist('doholdout','var'),       doholdout       = 0; end
if ~exist('donaivebayes','var'),    donaivebayes    = 1; end
if ~exist('doridge','var'),         doridge         = 0; end
if ~exist('dozscore','var'),        dozscore        = 1; end % ridge regression script does zscoring on its own.
if ~exist('dopermutelabels','var'), dopermutelabels = 0; end

timestep    = 'all';    %default: string 'all' to use all features, otherwise specify stepsize in sensors*samples.

nleaveout   = 2;        %if smaller than 1 ( but positive) will be interpreted as percentage      

folds      = 10;       % select number of folds for cross-validation

repetitions = 0;        % default no repetitions

%nfeatures   = 10000    % set how many features of gaussian naive bayes to
%use for testing
%% Load in & manipulate data
load(strcat(pathname));

[nobs,ncol]     = size(data);

if ischar(timestep), timestep = ncol; end

if dopermutelabels
    
    labels = labels(randperm(length(labels)));

end
%% Loop through repetitions, groups of features (sensors) & time windows & folds

acc = zeros(repetitions+1,length(timestep));

for nrep = 1:(1+repetitions)
        
        %%% CROSS-VALIDATION %%%
    
        if dofolds
        %stratified cross-validation 
            test_id = cvpartition(labels,'KFold',folds) 
        end
        
        if doholdout
        %random leave x out cross-validation
            test_id = cvpartition(nobs,'HoldOut',nleaveout) 
        end
        
        for ntime = timestep
            % select data according to features (i.e. time window)
            if length(ntime) == 1
                datap = data(:,1:ntime);
            elseif ntime+timestep < nobs
                datap = data(:,ntime:ntime+timestep);
            else 
                datap = data(:,ntime:end);
            end
            
            correct=0; incorrect=0;  rankAcc=[];  % counters for test accuracy
            
            for nfold = 1:folds
                % split trials into train & test
                data_test   = datap(test_id.test(nfold),:);
                data_train  = datap(test_id.training(nfold),:);
                label_test  = labels(test_id.test(nfold));
                label_train = labels(test_id.training(nfold));
                
                if dozscore
                   
                   [data_train, Mu, Sigma] = zscore(data_train);
                    data_test              = bsxfun(@rdivide, bsxfun(@minus, data_test, Mu),Sigma);
                
                end
                
                %%% TRAIN & PREDICT %%%
                
                if doridge
                    
                    [y_hat,beta_hat,lambda_hat]     = ridgeregression_sa(data_train,data_test,label_train,5,10);
                
                end
                
                if donaivebayes
                    
                    GNBmodel = nbayes_train(data_train, label_train, 1); 
        
                    % use top N features
                    %GNBmodel.sortedFeats   = sortGNBfeaturesByMuDifference(GNBmodel); 
                    %selectedFeats          = GNBmodel.sortedFeats(1:nFeatures);
                    
                    probs    = nbayes_apply(data_test, GNBmodel); 
        
                
                end
            
                %%% EVALUATE %%%
                if doridge && doholdout
                    
                    success = eval_euclideandistance(y_hat,y_test);
                    
                    if success == 0
                        fail = 1;
                    end
                    
                end
                
                if donaivebayes
                
                    [success,fail,rankA] = preps_decoding_eval_classprob(probs,label_test,GNBmodel,vocab);
       
                    rankAcc         = [rankAcc; rankA];

                end
                
                correct     = correct + success;
                incorrect   = incorrect + fail;
            end %end loop folds
            
            acc(nrep,size(ntime))=correct/(correct+incorrect);
            
        end %end loop time windows
end %end loop repetitions

