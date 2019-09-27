%In this script I inspect correlations between word2vec embedding vectors 

%load example dataset
load('/home/language/sopara/Projects/Prepositionalphrases/preps/Stimuli/preps_w2v.mat')

ind_noun = [];
ind_verb = [];
ind_adj  = [];
ind_det  = [];
ind_prep = [];
for i = 1:length(w2v.word)
    if any(strcmp({'NN'},w2v.pos(i)))
        ind_noun    = [ind_noun i];
    elseif any(strcmp({'VVFIN'},w2v.pos(i)))
        ind_verb    = [ind_verb i];
    elseif any(strcmp({'ADJA'},w2v.pos(i)))
        ind_adj    = [ind_adj i];
    elseif any(strcmp({'ART'},w2v.pos(i)))
        ind_det    = [ind_det i];
    elseif any(strcmp({'APPR'},w2v.pos(i)))
        ind_prep    = [ind_prep i];
    end
end

indices = [ind_noun ind_verb ind_adj ind_det ind_prep]';
feat_sort = w2v.feat(indices,:);
imagesc(corr(feat_sort'))

%dimensionality reduction
[coefs, s]  = pca(feat_sort);
low_data        = s*coefs(1:100,:)';


