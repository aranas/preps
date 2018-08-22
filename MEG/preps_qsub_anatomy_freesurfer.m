function preps_qsub_anatomy_freesurfer(subject,scriptname)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here


% Fressurfer script2
shell_script      = strcat('/home/language/sopara/Prepositionalphrases/preps/MEG/preps_anatomy_',scriptname,'.sh');
mri_dir           = '/project/3011210.01/anatomy';
subject_dir       = subject;

% streams_anatomy_freesurfer2.sh
command = [shell_script, ' ', mri_dir, ' ', subject_dir];

system(command);

end