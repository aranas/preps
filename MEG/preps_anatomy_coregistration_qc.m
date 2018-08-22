function preps_anatomy_coregistration_qc(save_dir, subject)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

headmodel   = fullfile(save_dir, subject, [subject '_headmodel.mat']);
sourcemodel = fullfile(save_dir, subject, [subject, '_sourcemodel.mat']);

load(headmodel)
load(sourcemodel)

figure; hold on;
ft_plot_vol(headmodel, 'facecolor', 'none'); alpha 0.5;
ft_plot_mesh(sourcemodel, 'facecolor', 'cortex', 'edgecolor', 'none'); camlight;

end

