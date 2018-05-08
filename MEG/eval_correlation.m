function [acc, dist, correct_pos, incorrect_pos] = eval_correlation(cfg,obj,testsize)

success = 0;
fail    = 0;
dist    = cell(length(cfg.nfolds),1);
incorrect_pos = {};
correct_pos    = {};

if testsize ~= size(obj.result{1},1)
    loop = size(obj.result{1},1)/testsize;
else
    loop = 1;
end
for f = 1:cfg.nfolds
    for l = 1:testsize:loop
        
        y_hat                = obj.result{f}(l:l+(testsize-1),:);
        y_test               = obj.design{f}(l:l+(testsize-1),:);
        
        [num,~]              = size(y_test);
        if mod(num,2)~= 0
            num = num-1;
        end
        y_test1              = y_test(1:(num / 2) , :)' ;
        y_test2              = y_test(num / 2 + 1:num , :)' ;
        y_test1              = y_test1(:) ;
        y_test2              = y_test2(:) ;
        
        y_hat1               = y_hat(1:num / 2 , :)' ;
        y_hat2               = y_hat(num/ 2 + 1:num , :)' ;
        y_hat1               = y_hat1(:) ;
        y_hat2               = y_hat2(:) ;
        
        %     correlations         = tril(corr([y_hat1 y_hat2 y_test1 y_test2]));
        %     correlations         = correlations - diag(diag(correlations));
        %     [m,n]                 = max(correlations(:));
        same                = corr(y_hat1,y_test1) + corr(y_hat2,y_test2);
        different           = corr(y_hat1,y_test2) + corr(y_hat2,y_test1);
        correct             = same > different;
        
        if correct
            success          = success + 1;
            correct_pos{f}   = obj.pos{f};
        else
            fail             = fail + 1;
            incorrect_pos{f}   = obj.pos{f};
        end
    end
end
acc = success / (success + fail);

correct_pos = correct_pos(~cellfun('isempty',correct_pos));
incorrect_pos = incorrect_pos(~cellfun('isempty',incorrect_pos));
end