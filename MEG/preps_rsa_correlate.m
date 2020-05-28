function out = preps_rsa_correlate(data,distancemeasure,model,correlation,perm_idx,normalize)

%%%%%%
% This function computes:
% 1. the dissimilarity matrix on the data looping over time points
% 2. the dissimilarity matrix for permuted trials
% 3. the correlation between all provided model RDMs and data RDM
%%%%%%
% input %%
% data is a vox_time data matrix
% distancemeasure used to compute data RDM (i.e. 'euclidean', 'correlation')
% model is a cell array struct with 
% model.rdm =  model RDMs
% model.name = corresponding condition names
% correlation used to correlate data and model RDMs (i.e. 'Spearman','Kendall')
% perm_idx contains the permuted trial ordering for each permutation.
n_trials = size(data,1);

for t = 1:size(data,3) %for feature (voxel or time)
    if normalize
        tmpdat = squareform(pdist(zscore(squeeze(data(:,:,t))')',distancemeasure));
    else
        tmpdat = squareform(pdist(squeeze(data(:,:,t)),distancemeasure));
    end
    for m = 1:length(model.rdm) %for each model
        out(m).name = model.name{m};
        out(m).rho(t) = corr(model.rdm{m}(tril(true(n_trials,n_trials),-1)),tmpdat(tril(true(n_trials,n_trials),-1)),'type',correlation,'rows', 'complete');
        
        for r = 1:size(perm_idx,1) %for each permutation
            tmpdat_r = tmpdat(perm_idx(r,:),perm_idx(r,:));
            out(m).rhorand(t,r) = corr(model.rdm{m}(tril(true(n_trials,n_trials),-1)),tmpdat_r(tril(true(n_trials,n_trials),-1)),'type',correlation,'rows', 'complete');
        end
    end
end