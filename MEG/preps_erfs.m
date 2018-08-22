% This script computes ERFs time-locked to the final word of each sentence
% It displayes the grand-average over all subjects
clear all;
clc;

trigger = {[111,114,121,124,211,214,221,224], %Determiner
    [112,115,122,125,212,215,222,225],        %Nouns
    [113,123,213,223],                        %Verbs
    [118,128,218,228],                        %Adjectives
    [116,126,216,226],                        %Preposition
    [219,229],                                %last word Noun attached
    [119,129],                                %last word Verb attached
    [30:39]};                                 %all words in filler sentences  
pos = {'Det','Noun','Verb','Adj','Prep','NA','VA','Fill'};
if ~exist('root_dir',    'var'), root_dir     = '/project/3011210.01/';  end
%% Load data & convert to planar gradients
filenames = '/project/3011210.01/MEG/sub*';
d = dir(filenames);
for subj = 1:length(d)
    load(strcat('/project/3011210.01/MEG/',d(subj).name))
    clear badcomp compds
  
    %reject noisy trials  
    cfg           = [];
    cfg.trials    =  ~ismember(data.trialinfo(:,3),noisy_trials(:,3));
    data          = ft_selectdata(cfg,data);
    
    %Baselining
    cfg           = [];
    cfg.demean    = 'yes';
    cfg.baselinewindow = [-0.1 0];
    data          = ft_preprocessing(cfg,data);
    
    %split conditions
    cfg           = [];
    cfg.trials    = ismember(data.trialinfo(:,1),[trigger{end-1}]);
    dataV         = ft_selectdata(cfg,data);
    cfg           = [];
    cfg.trials    = ismember(data.trialinfo(:,1),[trigger{end-2}]);
    dataN         = ft_selectdata(cfg,data);
    
    cfg           = [];
    cfg.removemean= 'yes';
    avgV{subj}          = ft_timelockanalysis(cfg, dataV);
    avgN{subj}          = ft_timelockanalysis(cfg, dataN);
    
    rawfiledir      = fullfile(sprintf('/project/3011210.01/raw/sub-%0.3d/ses-meg01/meg',subj));
    draw               = dir(rawfiledir);
    rawfile         = draw(3).name;
    sens{subj}      = ft_read_sens(fullfile(rawfiledir,rawfile),'senstype','meg');
end
%sens_avg=ft_average_sens([sens{:}]);FIXME: why is this not working?
for subj = 1:length(d)
    %Load headmodel and align to first subject 
    load(sprintf('/project/3011210.01/anatomy/sub-%0.3d/sub-%0.3d_headmodel.mat',subj,subj))
    
    cfg = [];
    cfg.headmodel   = headmodel;
    cfg.inwardshift = 2.5;
    cfg.template = sens{5};%sens_avg;
    cfg.feedback    = 'no';
    avgV_aligned{subj} = ft_megrealign(cfg, avgV{subj});
    avgN_aligned{subj} = ft_megrealign(cfg, avgN{subj});
    %Compute planar gradients
    cfg                 = [];
    cfg.feedback        = 'yes';
    cfg.method          = 'template';
    cfg.neighbours      = ft_prepare_neighbours(cfg, avgV_aligned{subj});
    
    cfg.planarmethod    = 'sincos';
    avgVplanar          = ft_megplanar(cfg, avgV_aligned{subj});
    avgNplanar          = ft_megplanar(cfg, avgN_aligned{subj});
    
    cfg = [];
    avgVplanarComb{subj} = ft_combineplanar(cfg,avgVplanar);
    avgNplanarComb{subj} = ft_combineplanar(cfg,avgNplanar);
    %calculate difference before or after GA?
end


cfg           = [];
gaV           = ft_timelockgrandaverage(cfg, avgVplanarComb{:});
gaN           = ft_timelockgrandaverage(cfg, avgNplanarComb{:});

save('/project/3011210.01/MEG/groupdata/group_erf_VAvsNA.m','gaV','gaN','sens');

%% Plotting
cfg = [];
cfg.xlim = [0.6 0.8];
% cfg.zlim = 'maxabs';
cfg.colorbar = 'yes';
cfg.layout = 'CTF275_helmet.mat';
ft_topoplotER(cfg,gaV,gaN)
colorbar;

%% Stats
cfg                   = [];
cfg.method            = 'distance';
cfg.neighbours        = ft_prepare_neighbours(cfg, avgV_aligned{5});
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

[stat] = ft_timelockstatistics(cfg, avgN_aligned{:}, avgV_aligned{:});