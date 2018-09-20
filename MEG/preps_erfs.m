% This script computes ERFs time-locked to the final word of each sentence
% It displayes the grand-average over all subjects

trigger = {[111,114,121,124,211,214,221,224], %Determiner
    [112,115,122,125,212,215,222,225],        %Nouns
    [113,123,213,223],                        %Verbs
    [118,128,218,228],                        %Adjectives
    [116,126,216,226],                        %Preposition
    [219,229],                                %last word Noun attached
    [119,129],                                %last word Verb attached
    [30:39]};                                 %all words in filler sentences
neg = {'Det','Noun','Verb','Adj','Prep','NA','VA','Fill'};
if ~exist('root_dir',     'var'), root_dir     = '/project/3011210.01/';  end
if ~exist('split_attach', 'var'), split_attach = 0;                       end
if ~exist('split_corpus', 'var'), split_corpus = 0;                       end
if ~exist('split_freq', 'var'),   split_freq = 1;                       end
%% Load data & convert to planar gradients
filenames = '/project/3011210.01/MEG/sub*dataclean_lp01.mat';
d = dir(filenames);
for subj = 1:length(d)
    load(strcat('/project/3011210.01/MEG/',d(subj).name))
    clear badcomp compds
    
    %reject noisy trials
%     load(fullfile(root_dir,'MEG',sprintf('sub-%0.3d_outlier',subj)));
%     ind_clean               = ~ismember(data.trialinfo(:,3),outlier_trials(:,end));
    cfg         = [];
    cfg.method  = 'summary';
    data        = ft_rejectvisual(cfg,data);
    %split conditions
    if split_attach
        cfg           = [];
        cfg.trials    = ismember(data.trialinfo(:,1),[trigger{end-1}]) & ind_clean;
        data1         = ft_selectdata(cfg,data);
        cfg           = [];
        cfg.trials    = ismember(data.trialinfo(:,1),[trigger{end-2}]) & ind_clean;
        data2         = ft_selectdata(cfg,data);
    elseif split_corpus
        cfg                   = [];
        cfg.trials            = ismember(data.trialinfo(:,1),[trigger{end-2} trigger{end-1}]) %& ind_clean;
        data                  = ft_selectdata(cfg,data);
        
        [midiff,data] = preps_getmi(data);     %better even:compute midiff at very beginning from stimulus file only, then reorder data according to 2nd column in trialinfo
        [val index] = sort(midiff);
        
        low = index(1:40);
        high = index(end-39:end);
        cfg           = [];
        cfg.trials    = false(size(data.trial));
        cfg.trials(low) = 1;
        data1         = ft_selectdata(cfg,data);
        cfg           = [];
        cfg.trials    = false(size(data.trial));
        cfg.trials(high) = 1;
        data2         = ft_selectdata(cfg,data);
        
    elseif split_freq
        cfg                   = [];
        cfg.trials            = ismember(data.trialinfo(:,1),[trigger{2:4} trigger{end-2} trigger{end-1}]) & ind_clean;
        data                  = ft_selectdata(cfg,data);
        
        load('/project/3011210.01/semanticP600/preps_stimuli.mat')  %load stimulus matrix to identify words
        %stimuli = stimuli([stimuli.condition]~=3);
        wikifreq = zeros(1,length(data.trial));
        for i = 1:length(data.trial)
            id = data.trialinfo(i,2);
            num = num2str(data.trialinfo(i,1));
            num = str2num(num(end));
            wikifreq(i)     = diag(stimuli(id).cooc(num,num));
        end

        [val index] = sort(wikifreq);
        
        low = index(1:round(length(wikifreq)/3));
        high = index(end-(round(length(wikifreq)/3-1)):end);
        cfg           = [];
        cfg.trials    = false(size(data.trial));
        cfg.trials(low) = 1;
        data1         = ft_selectdata(cfg,data);
        cfg           = [];
        cfg.trials    = false(size(data.trial));
        cfg.trials(high) = 1;
        data2         = ft_selectdata(cfg,data);
        
    end
    
    %Baselining
%     cfg           = [];
%     cfg.demean    = 'yes';
%     cfg.baselinewindow = [-0.2 0];
%     data1          = ft_preprocessing(cfg,data1);
%     data2          = ft_preprocessing(cfg,data2);
    
    cfg           = [];
    cfg.removemean= 'yes';
    cfg.vartrllength = 2;
    avg1{subj}          = ft_timelockanalysis(cfg, data1);
    avg2{subj}          = ft_timelockanalysis(cfg, data2);
    
    rawfiledir      = fullfile(sprintf('/project/3011210.01/raw/sub-%0.3d/ses-meg01/meg',subj));
    draw               = dir(rawfiledir);
    rawfile         = draw(3).name;
    sens{subj}      = ft_read_sens(fullfile(rawfiledir,rawfile),'senstype','meg');
end
%sens_avg=ft_average_sens([sens{:}]);FIXME: not working due to different
%sensors after MEG maintenance
for subj = 1:length(avg1)
    %Load headmodel and align to first subject
    load(sprintf('/project/3011210.01/anatomy/sub-%0.3d/sub-%0.3d_headmodel.mat',subj,subj))
    
    cfg = [];
    cfg.headmodel   = headmodel;
    cfg.inwardshift = 2.5;
    cfg.template = sens{1};%sens_avg;
    cfg.feedback    = 'no';
    avg1_aligned{subj} = ft_megrealign(cfg, avg1{subj});
    avg2_aligned{subj} = ft_megrealign(cfg, avg2{subj});
    %Compute planar gradients
    cfg                 = [];
    cfg.feedback        = 'yes';
    cfg.method          = 'distance';
    cfg.neighbours      = ft_prepare_neighbours(cfg, avg1_aligned{subj});
    
    cfg.planarmethod    = 'sincos';
    avg1planar          = ft_megplanar(cfg, avg1_aligned{subj});
    avg2planar          = ft_megplanar(cfg, avg2_aligned{subj});
    
    cfg = [];
    avg1planarComb{subj} = ft_combineplanar(cfg,avg1planar);
    avg2planarComb{subj} = ft_combineplanar(cfg,avg2planar);
    %calculate difference before or after GA?
    avg1planarComb{subj}.grad  = avg1_aligned{subj}.grad;  % add the gradiometer structure
    avg2planarComb{subj}.grad  = avg2_aligned{subj}.grad;
end


cfg           = [];
ga1           = ft_timelockgrandaverage(cfg, avg1planarComb{:});
ga2           = ft_timelockgrandaverage(cfg, avg2planarComb{:});

save('/project/3011210.01/MEG/groupdata/group_erf_lowfreqvshighfreq_allwords','gaV','gaN','sens','avg1planarComb','avg2planarComb');

%% Plotting
cfg = [];
cfg.xlim = [0.7 0.8];
%cfg.zlim = 'maxabs';
cfg.colorbar = 'yes';
cfg.layout = 'CTF275_helmet.mat';
figure()
ft_topoplotER(cfg,ga2,ga1)


%% Stats

cfg                 = [];
cfg.latency         = [-inf 1];
cfg.avgovertime     = 'no';
cfg.method          = 'distance';
cfg.neighbours      = ft_prepare_neighbours(cfg, avg1planarComb{1});
cfg.channel           = {'MEG'};
cfg.method            = 'montecarlo';
cfg.statistic         = 'depsamplesT';
cfg.correctm          = 'cluster';
cfg.clusteralpha      = 0.05;
cfg.clusterstatistic  = 'maxsum';
cfg.minnbchan         = 2;
cfg.tail              = 0;
cfg.clustertail       = 0;
cfg.alpha             = 0.025;
cfg.numrandomization  = 500;

subj                  = 10;
design                = zeros(2,2*subj);
for i = 1:subj
    design(1,i) = i;
end
for i = 1:subj
    design(1,subj+i) = i;
end
design(2,1:subj)        = 1;
design(2,subj+1:2*subj) = 2;

cfg.design            = design;
cfg.uvar              = 1;
cfg.ivar              = 2;

[stat] = ft_timelockstatistics(cfg, avg2planarComb{:}, avg1planarComb{:});

% Then take the difference of the averages using ft_math
cfg  = [];
cfg.operation = 'subtract';
cfg.parameter = 'avg';
highFvslowF = ft_math(cfg,ga2,ga1);

% Remember to do the same for negative clusters if you want them!

% First ensure the channels to have the same order in the average and in the statistical output.
% This might not be the case, because ft_math might shuffle the order  
[i1,i2] = match_str(highFvslowF.label, stat.label); 
unique(stat.negclusterslabelmat(stat.mask))

[~, firstcol] = find(stat.negclusterslabelmat==1,1);
[~, lastcol] = find(stat.negclusterslabelmat==1,1,'last');

for k = firstcol:4:lastcol%length(unique(stat.negclusterslabelmat(stat.mask)))
    figure();
    cfg = [];
    cfg.xlim =[stat.time(firstcol) stat.time(lastcol)];
    %cfg.zlim = [-1.0e-13 1.0e-13];   
    [rows,~] = find(stat.mask==1);
     cfg.highlight = 'on';
     cfg.highlightchannel = rows;
     cfg.comment = 'xlim';
     cfg.commentpos = 'title';
     cfg.layout = 'CTF275_helmet.mat';
     ft_topoplotER(cfg, highFvslowF);
end

cfg = [];
cfg.style     = 'blank';
cfg.layout    = 'CTF275_helmet.mat';
cfg.highlight = 'on';
cfg.highlightchannel = find(stat.mask);
cfg.comment   = 'no';
figure; ft_topoplotER(cfg, avg2planarComb{1})
title('Nonparametric: significant with cluster multiple comparison correction')