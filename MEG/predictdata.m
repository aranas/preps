function y_hat = predictdata(betas, mu, sigma, x)

x               = bsxfun(@rdivide, bsxfun(@minus, x, mu),sigma);
y_hat                = x * betas;
end