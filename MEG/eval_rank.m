function rankacc = eval_rank(cfg,obj)
% This script compares decoded vectors against all true vectors (by means
% of correlation) and computes the rank of the correct stimulus
% corresponding to each prediction. Ranks are then averaged across all
% decoded vectors and normalized to rank accuracy score (1- (rank-1/number of vectors in range -1))
% as described in Pereira et al. 2018

%combine predictions across folds
y_hat = cell2mat(obj.result)';
y_test = cell2mat(obj.design)';
%correlate matrix of predicted word embeddings with matrix of actual word
%embeddings
%de-mean both data and predicted data
C_1 = bsxfun(@minus, y_hat, mean(y_hat));
C_2 = bsxfun(@minus, y_test, mean(y_test));
%compute correlation between data & predicted data
R = corr(C_1,C_2);%sum(C_1 .* C_2) ./ (sqrt(sum(C_1 .^ 2)) .* sqrt(sum(C_2 .^ 2)));
%R has dimensions:npredicted X nactual
%rank
n = size(R,1);
rankR = zeros(size(R));
for vec = 1:n
   rankR(vec,:) = tiedrank(R(vec,:)); 
end
%average rank of correct stimulus for all items
avgrank = mean(diag(rankR));

%compute rank accuracy score
rankacc = 1-((avgrank-1)/(n-1));