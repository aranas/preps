function [acc, dist] = eval_euclideandistance(cfg,obj)

success = 0;
fail    = 0;
dist    = cell(length(cfg.nfolds),1);
for f = 1:cfg.nfolds
    
    y_hat                = obj.result{f};
    y_test               = obj.design{f};
    
    [num,~]              = size(y_test);
    y_test1              = y_test(1:num / 2 , :)' ;
    y_test2              = y_test(num / 2 + 1:end , :)' ;
    y_test1              = y_test1(:) ;
    y_test2              = y_test2(:) ;
    
    %evaluate prediction
    %euclidean distance between vectors belonging together should be smaller
    %than distance between vectors not belonging together.
    y_hat1               = y_hat(1:num / 2 , :)' ;
    y_hat2               = y_hat(num/ 2 + 1:end , :)' ;
    % number of trials is odd, therfore folds are odd, therefore test1
    % longer than test2 - problematic?
    y_hat1               = y_hat1(:) ;
    y_hat2               = y_hat2(:) ;
    same                 = (norm(y_hat1 - y_test1) + norm(y_hat2 - y_test2)) / 2;
    different            = (norm(y_hat1 - y_test2) + norm(y_hat2 - y_test1)) / 2;
    correct              = same < different;
    dist{f}              = [same different];
    
    if correct == 1
        success          = success + 1;
    else
        fail             = fail + 1;
    end
    acc = success / (success + fail);
    
end

end