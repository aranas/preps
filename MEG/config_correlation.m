%options for running preps_decoding inspired by discussion with Laura
%Gwilliams
%September 24th 2019

maincfg = [];
maincfg.mode = 'normal';
maincfg.classifier = 'ridgeregression_sa';
maincfg.seltrig = {'NN','VVFIN','ADJA'};
maincfg.dow2v = 1;
maincfg.w2vdim = 5;
maincfg.statfun = {'eval_correlation'};
maincfg.lambda = 1;
maincfg.toverlap = 0.8;
maincfg.repeats = 50;
%maincfg.dattype = 'lcmv';
%maincfg.parcel_indx = 1; %specify for each job
maincfg

suffix = '_eval2'