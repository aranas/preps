function success = eval_corr(y_hat,y_test,featv,num)
% which observations/words correspond to the test_set
success = 0;
d = size(y_hat);
 tmp = corr(y_test',featv');
 [maxval,testind] = max(tmp');
% correlate predicted vector with all possible vectors and choose top three
% picks.
 tmp = corr(y_hat',featv');
 [sorted,predictind] = sort(tmp,2,'descend');
 fprintf('highest correlation with any feature vectore was %.2f compared to 0.3376 by chance',max(sorted(:)))
 for i = 1:d(1)
 [row,col] = find(ismember(predictind(i,1:num),testind(i)));
 success = success + length(row);
 end
fprintf('%d correct matches within %d most correlated word embedding vectors',success,num);
fprintf('\n');
end