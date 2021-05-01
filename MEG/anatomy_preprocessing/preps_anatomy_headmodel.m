function preps_anatomy_headmodel(save_dir,subject)

headmodel_filename          = fullfile(save_dir, subject, [subject, '_headmodel' '.mat']);

mni_resliced_filename       = fullfile(save_dir, subject, [subject, '_mni_resliced' '.mgz']);
transform                   = fullfile(save_dir, subject, [subject, '_transform_vox2ctf.mat']);
load(transform);

mri                         = ft_read_mri(mni_resliced_filename);
mri.coordsys                = 'ctf';
mri.transform               = transform_vox2ctf;

cfg = [];
cfg.output = 'brain';
cfg.scalpthreshold = 0.9;
cfg.skullthreshold = 0.4;
cfg.brainthreshold = 0.4;
cfg.skullsmooth = 3; 
seg = ft_volumesegment(cfg, mri);

cfg = [];
cfg.method = 'projectmesh';
cfg.numvertices = 10000;
bnd = ft_prepare_mesh(cfg, seg);

cfg = [];
cfg.method = 'singleshell';
headmodel = ft_prepare_headmodel(cfg, bnd);

save(headmodel_filename, 'headmodel');

end

