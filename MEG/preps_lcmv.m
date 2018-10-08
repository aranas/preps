function [source_parc, filterlabel, source] = preps_lcmv(subj, data)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% computation of the covariance matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cfg              = [];
cfg.covariance   = 'yes';
cfg.vartrllength = 2;
cfg.channel      = 'MEG';
tlck = ft_timelockanalysis(cfg, data);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% preparation of the anatomical data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load the 2D sourcemodel and deal with the midline
load(fullfile('/project/3011210.01/anatomy',subj,strcat(subj,'_sourcemodel.mat')));

% define the medial wall parcel as outside. NOTE: this assumes
% the medial wall te have a value of 2
load atlas_subparc374_8k
exclude_label = match_str(atlas.parcellationlabel, {'L_???_01', 'L_MEDIAL.WALL_01', 'R_???_01', 'R_MEDIAL.WALL_01'});
sourcemodel.inside  = find(~ismember(atlas.parcellation,exclude_label));
sourcemodel.outside = find( ismember(atlas.parcellation,exclude_label));

% load the volume conduction model of the head
load(fullfile('/project/3011210.01/anatomy',subj,strcat(subj,'_headmodel.mat')));

% pre-compute the leadfields
cfg             = [];
cfg.grad        = ft_convert_units(tlck.grad,'m');
cfg.headmodel   = ft_convert_units(headmodel, 'm');
cfg.grid        = ft_convert_units(sourcemodel,'m');
cfg.channel     = 'MEG';
cfg.feedback    = 'textbar';
sourcemodel     = ft_prepare_leadfield(cfg, tlck);

cfg = [];
cfg.method = 'lcmv';
cfg.headmodel = headmodel;
cfg.grid      = sourcemodel;
cfg.lcmv.keepfilter = 'yes';
cfg.lcmv.fixedori   = 'yes';
cfg.lcmv.lambda     = '100%';
%cfg.lcmv.weightnorm = 'unitnoisegain';
source = ft_sourceanalysis(cfg, tlck);

F      = zeros(size(source.pos,1),numel(tlck.label));
F(source.inside,:) = cat(1,source.avg.filter{:});

% prepare the cfg for pca
cfg                       = [];
cfg.method                = 'pca';

tmp     = rmfield(data, {'elec' 'grad'});
selparc = setdiff(1:numel(atlas.parcellationlabel),exclude_label); % hard coded exclusion of midline and ???

source_parc.label = atlas.parcellationlabel(selparc);
source_parc.time  = tlck.time;
source_parc.F     = cell(numel(source_parc.label),1);
source_parc.avg   = zeros(numel(selparc),numel(source_parc.time));
source_parc.dimord = 'chan_time';

for k = 1:numel(selparc)
  tmpF = F(atlas.parcellation==selparc(k),:);
  tmp.trial = tmpF*data.trial;
  tmp.label = data.label(1:size(tmpF,1));
  tmpcomp   = ft_componentanalysis(cfg, tmp);

  source_parc.F{k}     = tmpcomp.unmixing*tmpF;
  source_parc.avg(k,:) = source_parc.F{k}(1,:)*tlck.avg;
end
filterlabel = tlck.label; % keep track of the channels that went into the spatial filters
