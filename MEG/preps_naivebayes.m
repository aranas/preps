function [beta_hat,probs,out] = preps_naivebayes(cfg,trainX,testX,trainY)

% [trainX, Mu, Sigma]     = zscore(trainX);
% testX                   = bsxfun(@rdivide, bsxfun(@minus, testX, Mu),Sigma);

 if ~isfield(cfg,'numfeat')
     cfg.numfeat = size(trainX,2);
 end

Mu = nanmean(trainX);
Sigma = nanstd(trainX);


idx = ~isnan(Mu);
if any(idx)
    trainX(:,idx) = bsxfun(@minus,trainX(:,idx),Mu(idx));
    testX(:,idx) = bsxfun(@minus,testX(:,idx),Mu(idx));
end

idx = ~isnan(Sigma) & ~(Sigma==0);
if any(idx)
    trainX(:,idx) = bsxfun(@rdivide,trainX(:,idx),Sigma(idx));
    testX(:,idx) = bsxfun(@rdivide,testX(:,idx),Sigma(idx));
end

beta_hat               = nbayes_train(trainX, trainY', cfg.poolsigma); %third argument yes/no pool sigma?

% use top N features
beta_hat.sortedFeats   = sortGNBfeaturesByMuDifference(beta_hat,1); %2nd argument useSigma to normalize MuDifference yes/no
selectedFeats          = beta_hat.sortedFeats(1:cfg.numfeat);

probs                  = nbayes_apply(testX, beta_hat,selectedFeats);
out                    = [Mu Sigma];
% [~,indx_predict]        = max(probs,[],2);
% 
% y_hat                   = indx_predict;