function preps_anatomy_dicom2mgz(root_dir,save_dir,subject)
%streams_anatomy_dicom2mgz takes the the subject info data structure (or subject string as 'sXX') 
%   
%   Picks up the dicom files, reslices the image and creates a .mgz file (spm coordsyst)

% select the last dicom file in subject's mri directory
dicom_dir  = strcat(root_dir,subject,'/ses-mri01/');


dicom_subdir = dir(dicom_dir);
dicom_subdir = dicom_subdir(end).name;
dicom_list = dir(fullfile(dicom_dir, dicom_subdir)); % choose the second folder with anatomical dicoms
dicom_file = fullfile(dicom_dir, dicom_subdir, dicom_list(end).name);

% read in the dicom files
mri   = ft_read_mri(dicom_file);

% filename for saving
mgz_filename = fullfile(save_dir, subject, [subject, '_mri' '.mgz']);

% save the images in the mgz format
cfg             = [];
cfg.filename    = mgz_filename;
cfg.filetype    = 'mgz';
cfg.parameter   = 'anatomy';
ft_volumewrite(cfg, mri);

end

