function [rho_mean, rho_var] = eval_correlation(cfg,obj)

%correlate prediction across trials with actual value for each dimension
%separately
rho = zeros(obj.folds,size(obj.result{1},2));
for fold = 1:obj.folds
    rho(fold,:) = diag(preps_corr(obj.result{fold},obj.design{fold}))';
end
rho_var = var(rho);
rho_mean= mean(rho);