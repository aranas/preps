function [beta_hat,y_hat,lambda_hat,lambdas,Mu,Sigma] = ridgeregression_sa(cfg,x_train,x_test,y_train)

cfg.nfolds = 10; 
if ~isfield(cfg, 'numlambdas'), cfg.numlambdas = 10;  end

%zscore
[x_train, Mu, Sigma] = zscore(x_train);

%pre-compute data variance
varx                 = x_train * x_train';
lambdas          = [0.000000001,0.0000000001,0.000000001,0.00000001,0.0000001,0.000001,1,5,8];
% train ridge regression with several lambda values using cross-validation
if exist('lambdas','var')
[R, lambdas]         = get_R_and_lambda(x_train,varx, y_train, cfg.nfolds, length(lambdas),lambdas); 
else
    [R, lambdas]         = get_R_and_lambda(x_train,varx,y_train,cfg.nfolds, cfg.numlambdas); %k =  folds for crossval; 10 lambda values
end
% force to use regularization:
R = R(:,2:end);
lambdas = lambdas(2:end);
% select optimal lambda values
[m,n]                = size(y_train);
r_hat                = NaN(n, 1);
lambda_hat           = NaN(n, 1);

for feat   = 1 : n
    
    [r_hat(feat), I] = max(R(feat, 1:end));
    lambda_hat(feat) = lambdas(I);
    
end

% Compute beta values with optimal Lambda
C                    = unique(lambda_hat);
beta_hat             = NaN(m, n);
I                    = eye(m);

for lambda  = 1 : length(C)
    
    beta_hat(:, C(lambda) == lambda_hat) = (varx + C(lambda) * I) \ y_train(:, C(lambda) == lambda_hat);

end
    beta_hat =  x_train' * beta_hat;
y_hat               = predictdata(beta_hat,Mu,Sigma,x_test);
end