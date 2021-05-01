function preps_anatomy_mgz2mni(anatomy_preproc_dir,subject)
% preps_anatomy_mgz2mni reads in the files created by
% preps_anatomy_dicom2mgz, and creates the transformation matrix for conversion to MNI space
% perfroms reslicing to 256x256x256 space. It saves the reliced volume and the transformation matrix. 

%% Initialize the variables


mgz_filename          = fullfile(anatomy_preproc_dir, subject, [subject '_mri.mgz']); 
resliced_filename     = fullfile(anatomy_preproc_dir, subject, [subject '_mni_resliced.mgz']);

% check if the .mgz file exists, if not create on spot
if ~exist(mgz_filename, 'file')
  
  fprintf('No .mgz found in %s\n', anatomy_preproc_dir);
  fprintf('Creating .mgz file via streams_anatomy_dicom2mgz\n');
  
  streams_anatomy_dicom2mgz(subject)

end

% load the mgz file
mri = ft_read_mri(mgz_filename);

%% Transform to MNI and reslice to freesurfer-friendly dimensions

% realign to MNI space
cfg                 = [];
cfg.coordsys        = 'spm';
cfg.parameter       = 'anatomy';
cfg.method          = 'interactive';
mri_mni             = ft_volumerealign(cfg, mri);

% reslice & save the transformation matrix to the anatomy_dir
cfg                 = [];
cfg.resolution      = 1;
cfg.dim             = [256 256 256];
mri_resliced        = ft_volumereslice(cfg, mri_mni);

% Save the resliced mni-transformed mri image
cfg                 = [];
cfg.filename        = resliced_filename;
cfg.filetype        = 'mgz';
cfg.parameter       = 'anatomy';
ft_volumewrite(cfg, mri_resliced)

% Save the transformation matrix
transform_vox2mni   = mri_resliced.transform;
filename_vox2mni    = fullfile(anatomy_preproc_dir, subject, [subject, '_transform_vox2mni']);
save(filename_vox2mni, 'transform_vox2mni');

end

