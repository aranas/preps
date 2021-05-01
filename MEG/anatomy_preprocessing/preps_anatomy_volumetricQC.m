function preps_anatomy_volumetricQC(anatomy_dir,subject)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

inp_dir         = fullfile(anatomy_dir, subject);

t1              = fullfile(inp_dir, 'mri', 'T1.mgz');
normalization2  = fullfile(inp_dir, 'mri', 'brain.mgz');
white_matter    = fullfile(inp_dir, 'mri', 'wm.mgz');
white_matter_old = fullfile(inp_dir, 'mri', 'wm_old.mgz');

% Show T1
mri = ft_read_mri(t1);
cfg = [];
cfg.interactive = 'yes';
ft_sourceplot(cfg, mri);
set(gcf, 'name', [subject ' ' 'T1'], 'numbertitle', 'off');

% Show skullstripped image
mri = ft_read_mri(normalization2);
cfg = [];
cfg.interactive = 'yes';
ft_sourceplot(cfg, mri);
set(gcf, 'name', [subject ' ' 'skull-stripped'], 'numbertitle', 'off');

% Show white matter image
mri = ft_read_mri(white_matter);
cfg = [];
cfg.interactive = 'yes';
ft_sourceplot(cfg, mri);
set(gcf, 'name', [subject ' ' 'white matter'], 'numbertitle', 'off');

if exist(white_matter_old)
  
  mri = ft_read_mri(white_matter_old);
  cfg = [];
  cfg.interactive = 'yes';
  ft_sourceplot(cfg, mri);
  set(gcf, 'name', [subject ' ' 'white matter old'], 'numbertitle', 'off');
  
end

end

