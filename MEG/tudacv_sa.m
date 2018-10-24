function [models, results, design] = tudacv_sa(X,Y,T,cfg,Gamma)
% Adapted from Vidaurre's toolbox https://github.com/OHBA-analysis/HMM-MAR/wiki/User-Guide
% changed to work with models other than regression


% Performs cross-validation of the TUDA model, which can be useful for
% example to compare different number or states or other parameters
% Cross-validation has a limitation for TUDA that is acknowledged in
%   Vidaurre et al (2017). ***TO COMPLETE***
%   (section 'Modelling between-trial temporal differences improves decoding
%   performance')
% NOTE: Specificying Nfeatures will lead to a circular assessment
% (overfitting of R2), so this is discouraged
%
% INPUT
% X: Brain data, (time by regions)
% Y: Stimulus, (time by q); q is no. of stimulus features
% T: Length of series
% options: structure with the training options - see documentation in
%                       https://github.com/OHBA-analysis/HMM-MAR/wiki
%  Apart from the options specified for tudatrain, specific to tudacv are:
%  - options.mode, referring to how to do cross-validation.
%                  . if 1, CV is done time point by time point, across CV folds,
%                   to produce a (time by 1) accuracy vector.
%                  . if 2, CV is performed across time points, so we would get
%                   a (time by time) accuracy "generalisation" matrix.
%  - options.CVmethod, cross-validation with TUDA has a fundamental issue:
%                  that the state time courses are estimated
%                  using both the data and the stimuli; then in the
%                  held-out data there is no way to know which state is
%                  active because, by definition, we do not use the
%                  stimulus (that is what we want to predict).
%                  This options establishes how to compute the model time
%                  courses in the held-out data
%                  . if 2, the state time courses are estimated using
%                  just data and linear regression, i.e. we try to predict
%                  the state time courses in held-out trials using the data
%                  . if 1, then it uses the state active at each time point is taken
%                  to be the one that is most active in training at each
%                  time point (losing between-trial variability)
%                  . if -1, then the state active at each time point is taken
%                  to be the one that is most active in held-out trials at each
%                  time point (thus also losing between-trial variability)
%                  . if -2, then the actual held-out state time courses are used.
%  - options.NCV, containing the number of cross-validation folds (default 10)
%  - options.lambda, regularisation penalty for the decoding
%  - options.lossfunc, loss function to compute the cross-validated error,
%                  default is 'quadratic', in which case R2 corresponds to
%                  explained variance. Other options are: 'absolute', 'huber'
%  - options.Nperm, if options.mode==2, it generate surrogates by permuting
%                   the state time courses across trials, such that the
%                   average remains the same but the trial-specific
%                   temporal features are lost
% Gamma: Precomputed decoding models time courses (optional)
%
% OUTPUT
% R2: cross-validated explained variance
% c: CV folds structure (c.training, c.test)
% pval: if options.Nperm > 1 and options.mode==2, this is the pvalue for
%       testing the state tiem courses: are these significantly meaningful,
%       above the average accuracy?
% surrogates: if options.Nperm > 1 and options.mode==2, these are the
%       surrogates created by permuting the state time courses across
%       trials
%
%
% Author: Diego Vidaurre, OHBA, University of Oxford (2017)
mvafun    = str2fun(cfg.mva.method{1});
N = length(T); q = size(Y,2); ttrial = T(1); p = size(X,2);
responses = permute(reshape(Y,[ttrial N q]),[2 3 1]);
responses = responses(:,:,1);
pval = [];

if isfield(cfg,'NCV')
    NCV = cfg.NCV; cfg = rmfield(cfg,'NCV');
else
    if isfield(cfg,'c')
        NCV = length(cfg.c.training);
    else
        NCV = 10;
    end
end
if isfield(cfg,'lambda')
    lambda = cfg.lambda; cfg = rmfield(cfg,'lambda');
else, lambda = 0.0001;
end
if isfield(cfg,'lossfunc')
    lossfunc = cfg.lossfunc; cfg = rmfield(cfg,'lossfunc');
else, lossfunc = 'quadratic';
end
if isfield(cfg,'Nfeatures') && cfg.Nfeatures>0 && cfg.Nfeatures<p
    warning('Specifying Nfeatures can lead to a biased calculation of CV-R2')
end
if isfield(cfg,'mode')
    mode = cfg.mode; cfg = rmfield(cfg,'mode');
else, mode = 1;
end
if isfield(cfg,'CVmethod')
    CVmethod = cfg.CVmethod; cfg = rmfield(cfg,'CVmethod');
else, CVmethod = 1;
end
if isfield(cfg,'Nperm')
    Nperm = cfg.Nperm; cfg = rmfield(cfg,'Nperm');
else, Nperm = 1;
end
if isfield(cfg,'verbose')
    verbose = cfg.verbose; cfg = rmfield(cfg,'verbose');
else, verbose = 1;
end
if nargin < 5, Gamma = []; end

if ~all(T==T(1)), error('All elements of T must be equal'); end

% Form CV folds; if response are categorical, then it's stratified
samples_per_value = length(Y) / length(unique(Y));
stratified = samples_per_value > 200;
if ~isfield(cfg,'c')
    if stratified
        %disp('Response is treated as categorical')
        tmp = zeros(length(responses),1);
        for j = 1:size(responses,2)
            rj = responses(:,j);
            uj = unique(rj);
            for jj = 1:length(uj)
                tmp(rj == uj(jj)) = tmp(rj == uj(jj)) + 100^(j-1) * jj;
            end
        end
        uy = unique(tmp);
        group = zeros(length(responses),1);
        for j = 1:length(uy)
            group(tmp == uy(j)) = j;
        end
        c2 = cvpartition(group,'KFold',NCV);
    else
        c2 = cvpartition(N,'KFold',NCV);
        %disp('Response is treated as continuous - no CV stratification')
    end
    c = struct();
    c.test = cell(NCV,1);
    c.training = cell(NCV,1);
    for icv = 1:NCV
        c.training{icv} = c2.training(icv);
        c.test{icv} = c2.test(icv);
    end; clear c2
else
    c = cfg.c; cfg = rmfield(cfg,'c');
end

% Preproc data and put in the right format
%[X,Y,T,cfg] = preproc4hmm(X,Y,T,cfg);
cfg = remove_options(cfg);
p = size(X,2); q = size(Y,2);
X = reshape(X,[ttrial N p]);
Y = reshape(Y,[ttrial N q]);
RidgePen = lambda * eye(p);

% Get state time courses
if isempty(Gamma)
    cfg.verbose = 0;
    [~,Gamma] = tudatrain(reshape(X,[ttrial*N p]),...
        reshape(Y,[ttrial*N q]),T,cfg);
end
K = size(Gamma,2);
Gamma = reshape(Gamma,[ttrial N K]);

% Estimate testing state time courses if necessary
if CVmethod > -2
    Gammapred = zeros(ttrial,N,K);
    for icv = 1:NCV
        Ntr = length(c.training{icv}); Nte = length(c.test{icv});
        Gammatrain = reshape(Gamma(:,c.training{icv},:),[ttrial Ntr K]);
        Gammatest = reshape(Gamma(:,c.test{icv},:),[ttrial Nte K]);
        switch CVmethod
            case 2
                Xtrain = permute(X(:,c.training{icv},:),[2 3 1]);
                Xtest = permute(X(:,c.test{icv},:),[2 3 1]);
                Gammatrain = permute(reshape(Gamma(:,c.training{icv},:),[ttrial Ntr K]),[2 3 1]);
                for t = 1:ttrial
                    B = (Xtrain(:,:,t)' * Xtrain(:,:,t) + 0.0001 * eye(p)) \ ...
                        Xtrain(:,:,t)' * Gammatrain(:,:,t);
                    pred = Xtest(:,:,t) * B;
                    pred = pred - repmat(min(min(pred,[],2), zeros(Nte,1)),1,K);
                    pred = pred ./ repmat(sum(pred,2),1,K);
                    Gammapred(t,c.test{icv},:) = pred;
                end
            case 1
                mGammatrain = squeeze(mean(Gammatrain,2));
                for t = 1:ttrial
                    [~,k] = max(mGammatrain(t,:));
                    Gammapred(t,c.test{icv},k) = 1;
                end
            case -1
                mGammatest = squeeze(mean(Gammatest,2));
                for t = 1:ttrial
                    [~,k] = max(mGammatest(t,:));
                    Gammapred(t,c.test{icv},k) = 1;
                end
        end
    end
else
    Gammapred = Gamma;
end

Ypred = zeros(ttrial,N,2,'single');%%FIXME: hard coded dimension

G = Gamma;
Gp = Gammapred;

for icv = 1:NCV
    Ntr = length(c.training{icv}); Nte = length(c.test{icv});
    Xtrain = reshape(X(:,c.training{icv},:),[ttrial*Ntr p]);
    ytrain = squeeze(Y(1,c.training{icv},:));
    Xtest = reshape(permute(X(:,c.test{icv},:),[2 1 3]),[Nte ttrial*p]);
    Gammatrain = reshape(G(:,c.training{icv},:),[ttrial*Ntr K]);
    Gammatest = reshape(Gp(:,c.test{icv},:),[ttrial*Nte K]);
    for k = 1:K
        sGamma = repmat(sqrt(Gammatrain(:,k)),1,p);
        Xtrain_k = Xtrain .* sGamma;
        Xtrain_k = reshape(permute(reshape(Xtrain_k,[ttrial Ntr p]),[2 1 3]),[Ntr ttrial*p]);
        %ytrain_k = ytrain .* sGamma(:,1:q);%% FIXME:valid to skip
        %this ste;??
        %Beta = (Xtrain_k' * Xtrain_k + RidgePen) \ (Xtrain_k' * ytrain_k);
        [model{k},result]  = mvafun(cfg,Xtrain_k,Xtest,ytrain');
        %sGamma = repmat(Gammatest(:,k),[1 q]);
        [maxval,maxind] = max(Gammatest,[],2);
        sGamma = (maxind==k);
        Ypred(:,c.test{icv},:) = Ypred(:,c.test{icv},:) + reshape(result .* sGamtuma, [ttrial Nte size(result,2)]);
    end
    results{icv} = reshape(Ypred(:,c.test{icv},:),[ttrial*Nte size(result,2)]);
    ytest = reshape(Y(:,c.test{icv},:),[ttrial*Nte q]);
    design{icv} = ytest;
    models{icv}.weights = model;
end
end


function options = remove_options(options)
% things that we don't want tudatrain to do again
if isfield(options,'filter'), options = rmfield(options,'filter'); end
if isfield(options,'Nfeatures'), options = rmfield(options,'Nfeatures'); end
if isfield(options,'detrend'), options = rmfield(options,'detrend'); end
if isfield(options,'onpower'), options = rmfield(options,'onpower'); end
if isfield(options,'standardise'), options = rmfield(options,'standardise'); end
if isfield(options,'embeddedlags'), options = rmfield(options,'embeddedlags'); end
if isfield(options,'pca'), options = rmfield(options,'pca'); end
options.add_noise = 0;
end


