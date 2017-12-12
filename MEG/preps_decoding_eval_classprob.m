


function [correct,incorrect,rankACC] = preps_decoding_eval_classprob(probs,label_test,GNBmodel,label_string)

rankACC             = rankAccuracy(probs,label_test,GNBmodel);

[~,indx_predict]    = max(probs,[],2);

label_hat           = GNBmodel.labelVocab(indx_predict);

correct = 0;
incorrect = 1;

fprintf('left: true labels - right: predicted labels ');
fprintf('\n');

for i = 1:size(label_hat,1)
    if label_hat(i) == label_test(i) 
        correct=correct+1;
        fprintf('correct! %d %s - %d %s',label_test(i),label_string{label_test(i)},label_hat(i),label_string{label_hat(i)});
    else
        incorrect=incorrect+1;
        fprintf('incorrect %d %s - %d %s',label_test(i),label_string{label_test(i)},label_hat(i),label_string{label_hat(i)});
    end
    
fprintf('\n');
       
end
        
                 