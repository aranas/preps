function [beta_hat,probs,out] = preps_naivebayes(cfg,trainX,testX,trainY)

% [trainX, Mu, Sigma]     = zscore(trainX);
% testX                   = bsxfun(@rdivide, bsxfun(@minus, testX, Mu),Sigma);

cfg.testonly     = ft_getopt(cfg, 'testonly','no');

if ~isfield(cfg,'numfeat')
    cfg.numfeat = size(trainX,2);
end

if strcmp('yes',cfg.testonly)
    if isfield(cfg,'model') && isfield(cfg,'out')
        beta_hat  = cfg.model;
        out       = cfg.out;
        testX     = bsxfun(@minus,testX,cfg.out.Mu);
        testX     = bsxfun(@rdivide,testX,Sigma);
    else
        warning('missing parameters, could be model or Mu/Sigma values')
    end
else
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
    out.Mu                 = Mu;
    out.Sigma              = Sigma;
end
% use top N features
[beta_hat.sortedFeats, vals]  = sortGNBfeaturesByMuDifference(beta_hat,1); %2nd argument useSigma to normalize MuDifference yes/no
out.Mudiff             = vals;
selectedFeats          = beta_hat.sortedFeats(1:cfg.numfeat);

probs                  = nbayes_apply(testX, beta_hat,selectedFeats);
% [~,indx_predict]        = max(probs,[],2);
%
% y_hat                   = indx_predict;