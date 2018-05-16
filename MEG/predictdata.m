function y_hat = predictdata(betas, mu, sigma, x,cfg)
if cfg.constant == 0
    x                = bsxfun(@rdivide, bsxfun(@minus, x, mu),sigma);
end
y_hat                = x * betas;
end