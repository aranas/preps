function success = eval_euclideandistance(y_hat,y_test)

[num,tmp]            = size(y_test);
y_test1              = y_test(1:num/2,:)' ;
y_test2              = y_test(num/2+1:end,:)' ;
y_test1              = y_test1(:) ;
y_test2              = y_test2(:) ;

%evaluate prediction
%euclidean distance between vectors belonging together should be smaller
%than distance between vectors not belonging together.
y_hat1               = y_hat(1:num/2,:)' ;
y_hat2               = y_hat(num/2+1:end,:)' ;
% number of trials is odd, therfore folds are odd, therefore test1
% longer than test2 - problematic?
y_hat1               = y_hat1(:) ;
y_hat2               = y_hat2(:) ;
same                 = (norm(y_hat1-y_test1)+ norm(y_hat2-y_test2))/2;
different            = (norm(y_hat1-y_test2)+ norm(y_hat2-y_test1))/2;
success              = same < different;
% fprintf('same = %.2f ', same);
% fprintf('different = %.2f', different);
% fprintf('distance test vectors = %.2f',norm(y_test1-y_test2))
% fprintf('\n');
end