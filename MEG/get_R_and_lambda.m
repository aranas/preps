function [R, lambda] = get_R_and_lambda(x_train,K, Y, k, n,lambda) % K = X_train * X_train' (covariance of training data), k = number of fold, n = number of lambda values

d       = size(Y);
%create new cross-validation set within training data for finding optimal lambda
Indices = sort(crossvalind('Kfold', d(1), k));
if ~exist('lambda','var')
     lambda  = get_lambda(K, n);
end

Y_hat   = NaN(d(1), d(2), n);

for index_1 = 1 : k % for each fold  
    
    
    Train = index_1 ~= Indices; % boolean vector
    Test  = index_1 == Indices; % boolean vector
    S     = sum(Train); %amount training examples
    N     = NaN(S, d(2), n);
    I     = eye(S);
    
    foo = K(Train, Train); % covariance of training set
    bar = Y(Train, :);
    
    parfor index_2 = 1 : n %for each lambda
        
        N(:, :, index_2) = (foo + lambda(index_2) * I) \ bar; % (XX' + Lambda*I) \ Y ---> 1st part of solution for beta weights, where: Betas = X_train'(X_trainX_train' + Lambda*I) \ Y_train

    end
    Y_hat(Test, :, :) = reshape(K(Test, Train) * reshape(N, S, d(2) * n), sum(Test), d(2), n); % ----> 2nd part of solution for Beta weights + prediction of Y_hat
                                                                                               % if Y = X * Betas or Y_test = X_test * Betas then after inserting Beta formula from above:
                                                                                               % Y = X_test*X_train'(X_train*X_train' + Lambda*I)\ Y_train, so that X_test*X_train' is K(Test,Train).
    % Y_hat contains predicted data for all observations as it is filled
    % for each fold
    
end

R = NaN(d(2), n);

for index = 1 : n
    %de-mean both data and predicted data
    C_1 = bsxfun(@minus, Y, mean(Y));
    C_2 = bsxfun(@minus, Y_hat(:, :, index), mean(Y_hat(:, :, index)));
    %compute correlation between data & predicted data
    R(:, index) = sum(C_1 .* C_2) ./ (sqrt(sum(C_1 .^ 2)) .* sqrt(sum(C_2 .^ 2)));
 
end

end

function lambda = get_lambda(K, n)

s      = svd(K);
s      = s(s > 0);
lambda = NaN(1, n);
L      = length(s);
df     = linspace(L, 1, n);
M      = mean(1 ./ s);

f       = @(df, lambda) df - sum(s ./ (s + lambda)     );
f_prime = @(    lambda)      sum(s ./ (s + lambda) .^ 2);

for index = 1 : n
    
    if index == 1
        
        lambda(index) = 0;
        
    else
        
        lambda(index) = lambda(index - 1);
        
    end
    
    lambda(index) = max(lambda(index), (L / df(index) - 1) / M);
    
    temp = f(df(index), lambda(index));
    
    tic;
    
    while abs(temp) > 1e-10
        
        lambda(index) = max(0, lambda(index) - temp / f_prime(lambda(index)));
        
        temp = f(df(index), lambda(index));
        
        if toc > 1
            
            warning('GET_LAMBDA did not converge.')
            
            break;
            
        end
        
    end
    
end

end