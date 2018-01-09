function [beta_hat,probs,Mu,Sigma] = preps_naivebayes(cfg,trainX,testX,trainY)

% [trainX, Mu, Sigma]     = zscore(trainX);
% testX                   = bsxfun(@rdivide, bsxfun(@minus, testX, Mu),Sigma);

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

beta_hat                = nbayes_train(trainX, trainY', 1);

% use top N features
%GNBmodel.sortedFeats   = sortGNBfeaturesByMuDifference(GNBmodel);
%selectedFeats          = GNBmodel.sortedFeats(1:nFeatures);

probs                   = nbayes_apply(testX, beta_hat);
% [~,indx_predict]        = max(probs,[],2);
% 
% y_hat                   = indx_predict;