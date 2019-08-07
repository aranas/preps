%This script uses the new MVPA toolbox for decoding to compare with my own code
subj = 'sub-003';
root_dir     = '/project/3011210.01/MEG/';
lcmvfile        = fullfile(root_dir,strcat(subj,'_preps_lcmv_parc.mat'));
channelfile     = fullfile(root_dir,sprintf('%s_dataclean%s.mat',subj));
artfctfile      = fullfile(root_dir,strcat(subj,'_muscle'));

seltrig = [219,229,119,129]; %[113 123 213 223 115 125 215 225];%
pos = {'NA','NA','VA','VA'};%{'VVFIN','VVFIN','VVFIN','VVFIN', 'NN','NN','NN','NN'};%pos = {'VVFIN', 'NN'};
if ~exist('pos','var') || isempty(seltrig)
[seltrig, pos] = preps_help_collecttrig(subj, seltrig);
end
load(channelfile,'data')

sel = ones(length(data.trialinfo),1);

load(artfctfile);
sel         = sel & ~ismember(data.trialinfo(:,3),noisy_trials(:,end));

sel             = sel & ismember(data.trialinfo(:,1),seltrig);
cfg             = [];
cfg.trials      = sel;
datasel         = ft_selectdata(cfg,data);
%fix:why does this step also slightly shift time axis?
% tn = size(datasel.time{1},2);
% for i = 1:length(datasel.trial)
%     datasel.time{i}(1:tn) = datasel.time{1}(1:tn);
% end

parcel_indx = 282;
num_comp = 5;
load(lcmvfile);
source_parc.filterlabel = filterlabel;
datasel = preps_sensor2parcel(datasel,source_parc,num_comp,parcel_indx);
clear source_parc source filterlabel

[~,loctrial]        = ismember(datasel.trialinfo(:,1),seltrig);
[upos, ulabel , labels]  = unique(pos(loctrial));

%Average over 100ms time snippets
twidth = 0.1;
toverlap = 0.8;
begtim = datasel.time{1}(1);
endtim = datasel.time{1}(end);

endtim = endtim - (toverlap*twidth);
time = linspace(begtim, endtim, round(abs(begtim-endtim) ./ ...
    (twidth - toverlap * twidth)) + 1);

stat = zeros(length(time),5);
for  t = 1:length(time)
rng('default'); 
cfg = [];
cfg.latency       = [time(t) time(t)+twidth];
cfg.avgovertime   = 'yes'; 
cfg.method = 'mvpa';
cfg.mvpa.classifier = 'lda';%goes fast for lda, but takes forever with svm
cfg.mvpa.metric = 'accuracy';
cfg.mvpa.k = 20;
cfg.mvpa.repeat = 50;
cfg.mvpa.stratify = 1;
%cfg.timextime = 'yes';
%cfg.latency = [0.27 0.37];
cfg.design = labels';
for rep = 1:5
out = ft_timelockstatistics(cfg,datasel); 
stat(t,rep)        = out.accuracy;
end
end
