%This script computes the RSA noise ceiling
%for corresponding RSA analysis on contentwords (pre final word)
clear all
all_seltrig     = {[118 128 218 228 112 122 212 222 125 215 225 215 113 123 213 223]};
datasuffix  = '_lp01';
suffix = '_overlap80';

toverlap = 0.8;
twidth = 0.1;
begtim = -0.2;
endtim = 0.8;
% For example parcel & each time slice
% parcel 117 {'L_22_B05_05'}
parcel_indx = 117;
jobid = {};
endtim = endtim - twidth;
time = linspace(begtim, endtim, round(abs(begtim-endtim) ./ ...
    (twidth - toverlap * twidth)) + 1);

datdist = cell(10,length(time));
for nsub = 1:10%For each subject load data and keep only distance matrix
    subj = sprintf('sub-%.3d',nsub);%map to full source space
    %% Load/select data & embeddings
    
    root_dir    = '/project/3011210.01/MEG';
    channelfile     = fullfile(root_dir,sprintf('%s_dataclean%s.mat',subj, datasuffix));
    lcmvfile        = fullfile(root_dir,strcat(subj,'_preps_lcmv_parc.mat'));
    load(channelfile,'data')
    
    for i = 1:length(all_seltrig)
        
        [seltrig{i}, pos{i}] = preps_help_collecttrig(subj, all_seltrig{i});
        sel = ones(length(data.trialinfo),1);
        sel             = sel & ismember(data.trialinfo(:,1),seltrig{i});
        cfg             = [];
        cfg.trials      = sel;
        datasel{i}         = ft_selectdata(cfg,data);
        
    end
    %sort according to story to keep order same across subjects
    [a,b] = sort(datasel{1}.trialinfo(:,2));
    datasel{1}.trialinfo = datasel{1}.trialinfo(b,:);
    datasel{1}.trial = datasel{1}.trial(b);
    datasel{1}.sampleinfo = datasel{1}.sampleinfo(b,:);
    datasel{1}.time = datasel{1}.time(b);
    
    clear data
    upos = uniqueStrCell(pos{end});
    suffix = [strcat(upos{:}) suffix];
    
    %convert to source data
    num_comp = 1;
    load(lcmvfile);
    source_parc.filterlabel = filterlabel;
    for i = 1:length(datasel)
        datasel{i} = preps_sensor2parcel(datasel{i},source_parc,num_comp,parcel_indx);
    end
    clear source_parc source filterlabel
    
    [brain_perword, ~, ~, ~, ~] = select_embeddings(datasel,pos,1);

    %% Create similarity matrices for data and Correlate with models
    %neural RDM per parcel (searchlight across time) & correlate
    for t = 1:length(time)
        btim = nearest(datasel{end}.time{1},time(t));
        etim = nearest(datasel{end}.time{1},time(t)+twidth);
        
        datdist{nsub,t} = squareform(pdist(zscore(squeeze(brain_perword(:,:,btim:etim))')','euclidean'));  
    end
end
% using kriegeskorte toolbox
ceiling_upperBound = zeros(size(datdist,2),1);
 for t = 1:size(datdist,2)
    [ceiling_upperBound(t), ceiling_lowerBound, bestFitRDM] = rsa.stat.ceilingAvgRDMcorr(cat(3,datdist{:,t}),'Spearman');
 end
%Split subject data (leave-one-out)
% subsel = 1:10;
% rho = nan(nsub,length(time));
% n_trials = size(datdist{1},1);
% for nsub = 1:size(datdist,1)
%     for t = 1:size(datdist,2)
%     lo_dat = datdist{nsub,t};
%     rest_dat = mean(cat(3,datdist{:,t}),3);
%     % For each split we compare left-out distance matrix to remaining
%     % distance matrices
%     
%     rho(nsub,t) = corr(rest_dat(tril(true(n_trials,n_trials),-1)),lo_dat(tril(true(n_trials,n_trials),-1)),'type','spearman','rows', 'complete');
%     end
% end
% % We average correlation coefficient over all splits
% upperCeil = mean(rho);
% %