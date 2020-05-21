% This script creates all figures for the thesis chapter on prepositional
% phrases based on scripts on github.com/aranas/preps, Commit XXXX

%% Results 1: Nouns vs Verbs 
%% send jobs per subject
clear all

maincfg             = [];
maincfg.mode        = 'normal';
maincfg.classifier  = 'preps_naivebayes';
maincfg.subj        = 'pilot-005';
maincfg.datasuffix  = '_lp01';
maincfg.twidth      = 0.1;
maincfg.toverlap    = 0.8; 
maincfg.seltrig     = [115,125,215,225,113,123,213,223];
maincfg.dattype     = 'sensor'; 
maincfg.dopca       = false;
maincfg.numfeat     = 'all';
maincfg.repeats     = 50;
suffix = '';


cd /project/3011210.01/tmp
for nsub = 1:10
    maincfg.subj = sprintf('sub-%.3d',nsub);
    qsubfeval('preps_execute_pipeline','preps_decoding',{'maincfg',maincfg},...
        'memreq',(1024^3)*10,'timreq',60*60*10,'batchid',sprintf('preps_timewindwos_100ms_%s',maincfg.subj))
end
maincfg.twidth      = 0.05;
for nsub = 1:10
    maincfg.subj = sprintf('sub-%.3d',nsub);
    qsubfeval('preps_execute_pipeline','preps_decoding',{'maincfg',maincfg},...
        'memreq',(1024^3)*10,'timreq',60*60*15,'batchid',sprintf('preps_timewindwos_50ms_%s',maincfg.subj))
end
maincfg.time      = 'auto';
for nsub = 1:10
    maincfg.subj = sprintf('sub-%.3d',nsub);
    qsubfeval('preps_execute_pipeline','preps_decoding',{'maincfg',maincfg},...
        'memreq',(1024^3)*10,'timreq',60*60*25,'batchid',sprintf('preps_timewindwos_allms_%s',maincfg.subj))
end

%when averaging over time
maincfg = rmfield(maincfg,'time');
maincfg.time_avg = true;
maincfg.twidth      = 0.1;

for nsub = 1:10
    maincfg.subj = sprintf('sub-%.3d',nsub);
    qsubfeval('preps_execute_pipeline','preps_decoding',{'maincfg',maincfg},...
        'memreq',(1024^3)*10,'timreq',60*60*10,'batchid',sprintf('preps_timewindwos_%s',maincfg.subj))
end
maincfg.twidth      = 0.05;
suffix      = strcat(suffix,'_50ms');
for nsub = 1:10
    maincfg.subj = sprintf('sub-%.3d',nsub);
    qsubfeval('preps_execute_pipeline','preps_decoding',{'maincfg',maincfg},{'suffix',suffix},...
        'memreq',(1024^3)*10,'timreq',60*60*15,'batchid',sprintf('preps_timewindwos_50ms_%s',maincfg.subj))
end

%pca & feature selection
clear all
maincfg             = [];
maincfg.mode        = 'normal';
maincfg.classifier  = 'preps_naivebayes';
maincfg.subj        = 'pilot-005';
maincfg.datasuffix  = '_lp01';
maincfg.twidth      = 0.1;
maincfg.toverlap    = 0.8; 
maincfg.seltrig     = [115,125,215,225,113,123,213,223];
maincfg.dattype     = 'sensor'; 
maincfg.dopca       = true;
maincfg.numfeat     = 'all';
maincfg.repeats     = 50;
suffix = '_pca';

maincfg.time_avg = false;
maincfg.save_full = true; 
for nsub = 1:10
    maincfg.subj = sprintf('sub-%.3d',nsub);
    qsubfeval('preps_execute_pipeline','preps_decoding',{'maincfg',maincfg},{'suffix',suffix},...
        'memreq',(1024^3)*10,'timreq',60*60*10,'batchid',sprintf('preps_timewindwos_%s',maincfg.subj))
end

maincfg.save_full = false; 
maincfg.twidth      = 0.05;
suffix      = strcat(suffix,'_50ms');
for nsub = 1:10
    maincfg.subj = sprintf('sub-%.3d',nsub);
    qsubfeval('preps_execute_pipeline','preps_decoding',{'maincfg',maincfg},{'suffix',suffix},...
        'memreq',(1024^3)*10,'timreq',60*60*15,'batchid',sprintf('preps_timewindwos_50ms_%s',maincfg.subj))
end

maincfg.time      = 'auto';
maincfg.twidth      = 0;
suffix = '_pca_allms';
for nsub = 1:10
    maincfg.subj = sprintf('sub-%.3d',nsub);
    qsubfeval('preps_execute_pipeline','preps_decoding',{'maincfg',maincfg},{'suffix',suffix},...
        'memreq',(1024^3)*10,'timreq',60*60*25,'batchid',sprintf('preps_timewindwos_allms_pca_%s',maincfg.subj))
end
maincfg = rmfield(maincfg,'time');

%when averaging over time
suffix = '_pca';
maincfg.time_avg = true;

for nsub = 1:10
    maincfg.subj = sprintf('sub-%.3d',nsub);
    qsubfeval('preps_execute_pipeline','preps_decoding',{'maincfg',maincfg},{'suffix',suffix},...
        'memreq',(1024^3)*10,'timreq',60*60*10,'batchid',sprintf('preps_timewindwos_%s',maincfg.subj))
end
maincfg.twidth      = 0.05;
suffix      = strcat(suffix,'_50ms');
for nsub = 1:10
    maincfg.subj = sprintf('sub-%.3d',nsub);
    qsubfeval('preps_execute_pipeline','preps_decoding',{'maincfg',maincfg},{'suffix',suffix},...
        'memreq',(1024^3)*10,'timreq',60*60*15,'batchid',sprintf('preps_timewindwos_50ms_%s',maincfg.subj))
end

%maincfg.pca = true;
maincfg.time_avg = false;
maincfg.numfeat = 150;
suffix = '_pca_featsel';

maincfg.twidth      = 0.1;
for nsub = 1:10
    maincfg.subj = sprintf('sub-%.3d',nsub);
    qsubfeval('preps_execute_pipeline','preps_decoding',{'maincfg',maincfg},{'suffix',suffix},...
        'memreq',(1024^3)*10,'timreq',60*60*10,'batchid',sprintf('preps_timewindwos_%s',maincfg.subj))
end
maincfg.twidth      = 0.05;
suffix      = strcat(suffix,'_50ms');
for nsub = 1:10
    maincfg.subj = sprintf('sub-%.3d',nsub);
    qsubfeval('preps_execute_pipeline','preps_decoding',{'maincfg',maincfg},{'suffix',suffix},...
        'memreq',(1024^3)*10,'timreq',60*60*15,'batchid',sprintf('preps_timewindwos_50ms_%s',maincfg.subj))
end

% different classifier methods
clear all
maincfg             = [];
maincfg.mode        = 'normal';
maincfg.classifier  = 'svm';
maincfg.subj        = 'pilot-005';
maincfg.datasuffix  = '_lp01';
maincfg.twidth      = 0.1;
maincfg.toverlap    = 0.8; 
maincfg.time_avg = false;
maincfg.seltrig     = [115,125,215,225,113,123,213,223];
maincfg.dattype     = 'sensor'; 
maincfg.dopca       = true;
maincfg.numfeat     = 'all';
maincfg.repeats     = 50;
suffix = '_pca';

for nsub = 1:10
    maincfg.subj = sprintf('sub-%.3d',nsub);
    qsubfeval('preps_execute_pipeline','preps_decoding',{'maincfg',maincfg},{'suffix',suffix},...
        'memreq',(1024^3)*10,'timreq',60*60*10,'batchid',sprintf('preps_svm_%s',maincfg.subj))
end

maincfg.classifier  = 'blogreg';
begtim = -0.2;
endtim = 0.8;
endtim = endtim - (maincfg.toverlap*maincfg.twidth);
time = linspace(begtim, endtim, round(abs(begtim-endtim) ./ ...
    (maincfg.twidth - maincfg.toverlap * maincfg.twidth)) + 1);
for st = 1:length(time)
    
    maincfg.time = [time(st)];

    suffix = sprintf('_pca_%i',st);
    for nsub = 1:10
        maincfg.subj = sprintf('sub-%.3d',nsub);
        qsubfeval('preps_execute_pipeline','preps_decoding',{'maincfg',maincfg},{'suffix',suffix},...
            'memreq',(1024^3)*10,'timreq',60*60*20,'batchid',sprintf('preps_blogreg_timeslice%i_%s',st,maincfg.subj))
    end
end

%collect blogreg per subject
rootdir = fullfile('/project','3011210.01','MEG','Classification');

for nsub = 1:10
    nsub
    subj = sprintf('sub-%.3d',nsub);
    allstat = {};
    allstatshuf = {};
    for st = 1:length(time)
        fname = dir(sprintf(fullfile(rootdir,'%s','sensor','blogreg_%s_lp01_20folds_*feats_NNVVFIN_pca_%i.mat'),subj,subj,st));
        load(fullfile(fname.folder,fname.name))
        allstat{st} = stat{1};
        allstatshuf = [allstatshuf; vertcat(statshuf(:))'];       
    end
    stat = allstat;
    statshuf = allstatshuf;
    cfgcv.time = time;
    save(sprintf(fullfile(rootdir,'%s','sensor','blogreg_%s_lp01_20folds_1860feat_NNVVFIN_pca_all.mat'),subj,subj),'stat','statshuf','cfgcv');
end


%word position control
clear all

maincfg             = [];
maincfg.mode        = 'normal';
maincfg.classifier  = 'preps_naivebayes';
maincfg.subj        = 'pilot-005';
maincfg.datasuffix  = '_lp01';
maincfg.twidth      = 0.1;
maincfg.toverlap    = 0.8; 
maincfg.seltrig     = [33,35];
maincfg.dattype     = 'sensor'; 
maincfg.dopca       = true;
maincfg.numfeat     = 'all';
maincfg.repeats     = 50;
suffix = '_pca_filler';

for nsub = 1:10
    maincfg.subj = sprintf('sub-%.3d',nsub);
     %overwrite labeling by providing pos variable (pretending every 3rd and 5th word in fillers is actually NN and VVFIN
    qsubfeval('preps_execute_pipeline','preps_decoding',{'maincfg',maincfg},{'suffix',suffix},{'pos',{'NN','VVFIN'}},...
        'memreq',(1024^3)*10,'timreq',60*60*15,'batchid',sprintf('preps_timewindwos_poscontrol_%s',maincfg.subj))
end

%generalise to final word
clear all
maincfg             = [];
maincfg.mode        = 'general';
maincfg.classifier  = 'preps_naivebayes';
maincfg.subj        = 'pilot-005';
maincfg.datasuffix  = '_lp01';
maincfg.twidth      = 0.1;
maincfg.toverlap    = 0.8;
maincfg.time        = [0 2.1];
maincfg.seltrig     = [115,125,215,225,113,123,213,223];
maincfg.dattype     = 'sensor'; 
maincfg.dopca       = true;
maincfg.numfeat     = 'all';
maincfg.repeats     = 50;
maincfg.save_full   = false;
maincfg.testtrig    = [119,129,219,229];
maincfg.testpos     = {'VVFIN','VVFIN','NN','NN'};


begtim = 0;
endtim = 0.8;
endtim = endtim - (maincfg.toverlap*maincfg.twidth);
time = linspace(begtim, endtim, round(abs(begtim-endtim) ./ ...
    (maincfg.twidth - maincfg.toverlap * maincfg.twidth)) + 1);

for st = 1:length(time)-1
    maincfg.trainwindow = [time(st) time(st)+maincfg.twidth];
    suffix = sprintf('_pca_%i',st);
    for nsub = 1:10
        maincfg.subj = sprintf('sub-%.3d',nsub);
        qsubfeval('preps_execute_pipeline','preps_decoding',{'maincfg',maincfg},{'suffix',suffix},...
        'memreq',(1024^3)*10,'timreq',60*60*10,'batchid',sprintf('preps_general_%s_t%i',maincfg.subj,st))
    end
end

%collect generalize
rootdir = fullfile('/project','3011210.01','MEG','Classification');

for nsub = 1:10
    subj = sprintf('sub-%.3d',nsub);
    for st = 1:length(time)-1
        fname = dir(sprintf(fullfile(rootdir,'%s','sensor','nbayes_%s_lp01_20folds*__pca_%i_general.mat'),subj,subj,st));
        load(fullfile(fname.folder,fname.name))
        stat = struct2cell(cell2mat(stat));
        statshuf = struct2cell(cell2mat(statshuf));
        acc(st,:,nsub) = cell2mat(squeeze(stat(1,:,:)));
        accshuf(st,:,:,nsub) = cell2mat(squeeze(statshuf(1,:,:)));
    end
end
cfgcv.trainwind = time(1:end-1);
save(fullfile(rootdir,'groupresults','NNVVFIN_generalize.mat'),'cfgcv','acc','accshuf');

% source space
clear all

maincfg             = [];
maincfg.mode        = 'normal';
maincfg.classifier  = 'preps_naivebayes';
maincfg.subj        = 'pilot-005';
maincfg.datasuffix  = '_lp01';
maincfg.twidth      = 0.1;
maincfg.toverlap    = 0.8; 
maincfg.seltrig     = [115,125,215,225,113,123,213,223];
maincfg.dattype     = 'lcmv'; 
maincfg.dopca       = false;
maincfg.numfeat     = 'all';
maincfg.repeats     = 50;
maincfg.save_full   = true;
suffix = '_new';

%send off jobs for computing parcel-wise results
for nsub = 1:10
    maincfg.subj = sprintf('sub-%.3d',nsub);
    for parcid = 1:370
        maincfg.parcel_indx = parcid;
        qsubfeval('preps_execute_pipeline','preps_decoding',{'maincfg',maincfg},{'suffix',suffix},...
        'memreq',(1024^3)*5,'timreq',60*60*5,'batchid',sprintf('preps_timewindwos_%s_lcmv_%i',maincfg.subj,maincfg.parcel_indx))
    end
    pause((60*60*2))
end

%combine parcel-wise results into one file
rootdir = fullfile('/project','3011210.01','MEG','Classification');

for nsub = 8:10
    subj = sprintf('sub-%.3d',nsub);
    acc = zeros(47,370);
    accshuf = zeros(47,50,370);
    for i = 1:370 %there should be 370 parcel per subject
        fname = dir(sprintf(fullfile(rootdir,'%s','lcmv','searchlight','*_parcel%.3d_new.mat'),subj,i));
        load(fullfile(fname.folder,fname.name))
        stat = struct2cell(cell2mat(stat));
        statshuf = struct2cell(cell2mat(statshuf));
        acc(:,i) = cell2mat(squeeze(stat(1,:,:)));
        accshuf(:,:,i) = cell2mat(squeeze(statshuf(1,:,:)));
        i
    end  
    save(fullfile(fname.folder,sprintf('%s_acc_allparcels_new.mat',subj)),'cfgcv','acc','accshuf')
end

%NA vs VA

clear all
maincfg             = [];
maincfg.mode        = 'normal';
maincfg.classifier  = 'preps_naivebayes';
maincfg.subj        = 'pilot-005';
maincfg.datasuffix  = '_lp01';
maincfg.time        = [-0.2 2.1];
maincfg.twidth      = 0.1;
maincfg.toverlap    = 0.8; 
maincfg.seltrig     = {'NA', 'VA'};
maincfg.dattype     = 'sensor'; 
maincfg.dopca       = true;
maincfg.numfeat     = 'all';
maincfg.repeats     = 50;
suffix = '_pca';

maincfg.time_avg = false;
maincfg.save_full = false; 
for nsub = 1:10
    maincfg.subj = sprintf('sub-%.3d',nsub);
    qsubfeval('preps_execute_pipeline','preps_decoding',{'maincfg',maincfg},{'suffix',suffix},...
        'memreq',(1024^3)*10,'timreq',60*60*10,'batchid',sprintf('preps_NAVA_%s',maincfg.subj))
end

maincfg.doposttest = 1;
for nsub = 1:10
    maincfg.subj = sprintf('sub-%.3d',nsub);
    qsubfeval('preps_execute_pipeline','preps_decoding',{'maincfg',maincfg},{'suffix',suffix},...
        'memreq',(1024^3)*10,'timreq',60*60*10,'batchid',sprintf('preps_NAVA_%s',maincfg.subj))
end

% maincfg.classifier = 'blogreg';
% begtim = -0.2;
% endtim = 2.1;
% maincfg.toverlap = 0.5;
% endtim = endtim - (maincfg.toverlap*maincfg.twidth);
% time = linspace(begtim, endtim, round(abs(begtim-endtim) ./ ...
%     (maincfg.twidth - maincfg.toverlap * maincfg.twidth)) + 1);
% for st = 1:length(time)
%     
%     maincfg.time = [time(st)];
% 
%     suffix = sprintf('_pca_%i',st);
%     for nsub = 1:10
%         maincfg.subj = sprintf('sub-%.3d',nsub);
%         qsubfeval('preps_execute_pipeline','preps_decoding',{'maincfg',maincfg},{'suffix',suffix},...
%             'memreq',(1024^3)*10,'timreq',60*60*20,'batchid',sprintf('preps_blogreg_NAVA_t%i_%s',st,maincfg.subj))
%     end
% end

%collect blogreg per subject
% rootdir = fullfile('/project','3011210.01','MEG','Classification');
% 
% for nsub = 1:10
%     nsub
%     subj = sprintf('sub-%.3d',nsub);
%     allstat = {};
%     allstatshuf = {};
%     for st = 1:length(time)
%         fname = dir(sprintf(fullfile(rootdir,'%s','sensor','blogreg_%s_lp01_20folds_*feats_NAVA_pca_%i.mat'),subj,subj,st));
%         load(fullfile(fname.folder,fname.name))
%         allstat{st} = stat{1};
%         allstatshuf = [allstatshuf; vertcat(statshuf(:))'];       
%     end
%     stat = allstat;
%     statshuf = allstatshuf;
%     cfgcv.time = time;
%     save(sprintf(fullfile(rootdir,'%s','sensor','blogreg_%s_lp01_20folds_1860feat_NAVA_pca_all.mat'),subj,subj),'stat','statshuf','cfgcv');
% end


% Send jobs for w2v ridge regression
clear all
maincfg = [];
maincfg.mode = 'normal';
maincfg.classifier = 'ridgeregression_sa';
maincfg.statistic = {'eval_correlation'};
maincfg.dow2v = 1;
maincfg.seltrig = [115,125,215,225,113,123,213,223,118,128,218,228,119,129,219,229];
maincfg.type = 'nfold';
maincfg.lambdaeval= 'mse';
maincfg.lambdas = [10 10*exp(2) 10*exp(4) 10*exp(6)]';
maincfg.resample = 0;
maincfg.twidth = 0.1;
maincfg.toverlap = 0.5;
maincfg.repeats = 50;

begtim = -0.2;
endtim = 0.8;
endtim = endtim - (maincfg.toverlap*maincfg.twidth);
time = linspace(begtim, endtim, round(abs(begtim-endtim) ./ ...
    (maincfg.twidth - maincfg.toverlap * maincfg.twidth)) + 1);

for st = 1:length(time)-1
    maincfg.trainwindow = [time(st) time(st)+maincfg.twidth];
    suffix = sprintf('_pca_%i',st);
    for nsub = 1:10
        maincfg.subj = sprintf('sub-%.3d',nsub);
        qsubfeval('preps_execute_pipeline','preps_decoding',{'maincfg',maincfg},{'suffix',suffix},...
        'memreq',(1024^3)*10,'timreq',60*60*10,'batchid',sprintf('preps_NAVAw2v_%s_t%i',maincfg.subj,st))
    end
end


%% collect results for different parameter settings
clear all
rootdir = fullfile('/project','3011210.01','MEG','Classification');
savedir = fullfile('/project','3011210.01','MEG','figures','final');

cmap = brewermap(6,'Pastel1'); %for plotting lines
cmap_sig = brewermap(4,'OrRd'); %for plotting significance

fnames = {'nbayes_%s_lp01_20folds_27*feats_NNVVFIN.mat'... %1) all time points, timeXsensor
    'nbayes_%s_lp01_20folds_35*feats_NNVVFIN.mat'...       %2) 50ms window, timeXsensors
    'nbayes_%s_lp01_20folds_*feats_NNVVFIN_50ms_avg.mat'...%3) 50ms window, average
    'nbayes_%s_lp01_20folds_67*feats_NNVVFIN.mat'...       %4) 100ms window, timeXsensor
    'nbayes_%s_lp01_20folds_27*feats_NNVVFIN_avg.mat'...   %5) 100ms window, average
    'nbayes_%s_lp01_20folds_1500feats_NNVVFIN_pca.mat'...  %6) 100ms, timeXsensor, pca
        'nbayes_%s_lp01_20folds_150feats_NNVVFIN_pca_featsel.mat'...  %7) 100ms, timeXsensor, pca, featsel
    'nbayes_%s_lp01_20folds_60feats_NNVVFIN_pca_avg.mat'...           %8) 100ms, average, pca, 
    'nbayes_%s_lp01_20folds_780feats_NNVVFIN_pca_50ms.mat'...         %9) 50ms, timeXsensor, pca,
    'nbayes_%s_lp01_20folds_150feats_NNVVFIN_pca_featsel_50ms.mat'... %10) 50ms, timeXsensor, pca, featsel
    'nbayes_%s_lp01_20folds_60feats_NNVVFIN_pca_50ms_avg.mat'...      %11) 50ms, average, pca
    'nbayes_%s_lp01_20folds_60feats_NNVVFIN_pca_allms.mat'...         %12) allT, timexsensor,pca
    'nbayes_%s_lp01_20folds_1500feats_NNVVFIN_pca_filler.mat'...      %13) Nouns vs Verbs word position control
    'nbayes_%s_lp01_20folds_1860feats_NAVA_pca.mat'...                %14) Noun-attached vs Verb-attached 
    'nbayes_%s_lp01_20folds_1860feats_NAVA_pca_posttest.mat'...       %15) NA vs VA after relabeling according to individual posttest
    'svm_%s_lp01_20folds_1500feats_NNVVFIN_pca.mat'...                 %16)Noun vs Verbs with SVM classifier
        'blogreg_%s_lp01_20folds_1860feat_NNVVFIN_pca_all.mat'};

for nsub = 1:10
    subj = sprintf('sub-%.3d',nsub);
    %the exact number of features in file name can vary as some subjects
    %had fewer sensors recorded (due to broken MEG sensors)
    for i = 1:length(fnames)
        filename = dir(sprintf(fullfile(rootdir,'%s','sensor',fnames{i}),subj,subj)); 
        load(fullfile(filename.folder,filename.name))
        stat = struct2cell(cell2mat(stat));
        statshuf = struct2cell(cell2mat(statshuf));
        cfg{i} = cfgcv;
        acc{i}(:,nsub) = cell2mat(squeeze(stat(1,:,:)));
        accshuf{i}(:,:,nsub) = cell2mat(squeeze(statshuf(1,:,:)));
        time{i} = cfgcv.time;
    end
end

%% Plot1: different feature transformations on 100 ms moving time window, while concatenating time X Sensors:
plotcfg = {[4,6,7]};
taxis = round(time{plotcfg{1}(1)}*1000);

%plot permutation distribution
maccshuf = mean(reshape(accshuf{plotcfg{1}(1)},length(taxis),[]),2);
accshufbounds = std(reshape(accshuf{plotcfg{1}(1)},length(taxis),[]),[],2);
figure('units','normalized','outerposition',[0 0 1 1]); hold on;
hl = [];
hl(1)= boundedline(taxis, maccshuf, accshufbounds, ...
    'alpha','cmap',cmap(end,:));

%plot data
for l = 1:length(plotcfg{1})
    datid = plotcfg{1}(l);
    
    macc = mean(acc{datid},2);
    accbounds = std(acc{datid},[],2);
    
    hl(1+l) = boundedline(taxis, macc, accbounds, ...
        'alpha','cmap',cmap(l,:));
end

%add legend & title
title('Accuracy for classification of nouns vs. verbs')
xlabel(sprintf('time in ms (center of %d ms sliding window)',cfg{plotcfg{1}(l)}.twidth*1000))
ax = gca;
L = [];
L(1) = plot(nan,nan,'.','MarkerEdgeColor',cmap(1,:),'MarkerSize',40,'Parent',ax);
L(2) = plot(nan,nan,'.','MarkerEdgeColor',cmap(2,:),'MarkerSize',40,'Parent',ax);
L(3) = plot(nan,nan,'.','MarkerEdgeColor',cmap(3,:),'MarkerSize',40,'Parent',ax);
L(4) = plot(nan,nan,'.','MarkerEdgeColor',cmap(end,:),'MarkerSize',40,'Parent',ax);
hl(end+1) = legend(L,{'raw features','pca','pca + feature selection','permuted class labels'},'Location','best');

%format
ylim([0.4 0.8])
xlim([-200 700])
for h = 1:length(hl)
    set(hl(h),'Linewidth',5)
end
set(hl(end),'FontSize',25)
set(ax,'FontSize',25)

print(gcf,fullfile(savedir,'plot1.eps'),'-depsc','-painters')

%% Plot2: different treatment of time dimension for pca-transformed data:
plotcfg = {[6,8],[9,11],[12]};
taxis = round(time{plotcfg{1}(1)}*1000);

%plot permutation distribution
maccshuf = mean(reshape(accshuf{plotcfg{1}(1)},length(taxis),[]),2);
accshufbounds = std(reshape(accshuf{plotcfg{1}(1)},length(taxis),[]),[],2);
figure('units','normalized','outerposition',[0 0 1 1]); hold on;
hl = [];ax = {};
hl(1)= boundedline(taxis, maccshuf, accshufbounds, ...
    'alpha','cmap',cmap(end,:));

%plot data
ax{1} = gca;
ax1_Pos = ax{1}.Position;
ccnt = 1;
for p = 1:length(plotcfg)
    taxis = round(time{plotcfg{p}(1)}*1000);
    ax{1+p} = axes('Position',ax1_Pos,'color','none');   
for l = 1:length(plotcfg{p})
    datid = plotcfg{p}(l);
    
    macc = mean(acc{datid},2);
    accbounds = std(acc{datid},[],2); % use nan(size(macc)) to not plot bounds for better visibility
    
    hl(1+ccnt) = boundedline(taxis, macc, accbounds, ax{1+p}, ...
        'alpha','cmap',cmap(ccnt,:));
    set(ax{1+p},'Visible','off') 
    ccnt = ccnt+1;
end
end

%add legend & title
title('Accuracy for classification of nouns vs. verbs','Parent',ax{1})
xlabel('time in ms (center of sliding window)','Parent',ax{1})
L = [];
L(1) = plot(nan,nan,'.','MarkerEdgeColor',cmap(1,:),'MarkerSize',40,'Parent',ax{1});
L(2) = plot(nan,nan,'.','MarkerEdgeColor',cmap(2,:),'MarkerSize',40,'Parent',ax{1});
L(3) = plot(nan,nan,'.','MarkerEdgeColor',cmap(3,:),'MarkerSize',40,'Parent',ax{1});
L(4) = plot(nan,nan,'.','MarkerEdgeColor',cmap(4,:),'MarkerSize',40,'Parent',ax{1});
L(5) = plot(nan,nan,'.','MarkerEdgeColor',cmap(5,:),'MarkerSize',40,'Parent',ax{1});
L(6) = plot(nan,nan,'.','MarkerEdgeColor',cmap(end,:),'MarkerSize',40,'Parent',ax{1});
hl(end+1) = legend(L,{'100ms window & timeXsensors','100ms window & average','50ms window & timeXsensors','50ms window & average','single time points','permuted class labels'});

%format
for h = 1:length(hl)
    set(hl(h),'Linewidth',5)
end
for a = 1:length(ax)
    ylim(ax{a},[0.45 0.75])
    xlim(ax{a},[-200 700]) 
end
set(hl(end),'FontSize',25)
set(ax{1},'FontSize',25)

print(gcf,fullfile(savedir,'plot2.eps'),'-depsc','-painters')


%% Plot3
plotcfg = {[6,17]};
taxis = round(time{plotcfg{1}(1)}*1000);
%stat

dat.dimord = 'chan_time';
dat.label = {'all'};
dat.time = time{plotcfg{1}(1)};
for nsub = 1:10
    acc_NV{nsub} = dat;
    acc_NV{nsub}.avg = acc{plotcfg{1}(1)}(:,nsub)';
    acc_pos{nsub} = dat;
    acc_pos{nsub}.avg = acc{plotcfg{1}(2)}(:,nsub)';  
end

cfgs = [];
cfgs.method = 'montecarlo';
cfgs.neighbours = [];
cfgs.statistic = 'ft_statfun_depsamplesT';
cfgs.correctm = 'cluster';
cfgs.correcttail = 'prob';
cfgs.clusteralpha = 0.05;
cfgs.clusterstatistic = 'maxsum';
cfgs.alpha = 0.025;
cfgs.numrandomization = 'all';
cfgs.design(1,1:2*nsub) = [ones(1,nsub) 2*ones(1,nsub)];
cfgs.design(2,1:2*nsub) = [1:nsub 1:nsub];
cfgs.ivar = 1;
cfgs.uvar = 2;
stat = ft_timelockstatistics(cfgs,acc_NV{:},acc_pos{:});

%plot permutation distribution
maccshuf = mean(reshape(accshuf{plotcfg{1}(1)},length(taxis),[]),2);
accshufbounds = std(reshape(accshuf{plotcfg{1}(1)},length(taxis),[]),[],2);
figure('units','normalized','outerposition',[0 0 1 1]); hold on;
hl = [];
hl(1)= boundedline(taxis, maccshuf, accshufbounds, ...
    'alpha','cmap',cmap(end,:));

%plot data
for l = 1:length(plotcfg{1})
    datid = plotcfg{1}(l);
    
    macc = mean(acc{datid},2);
    accbounds = std(acc{datid},[],2);
    
    hl(1+l) = boundedline(taxis, macc, accbounds, ...
        'alpha','cmap',cmap(l,:));
end
%Plot significance indicator
idx_sig = double(stat.mask);
idx_sig(idx_sig==0) = NaN;
plot(taxis,(0.5-0.06)*idx_sig,'-k','LineWidth',5)

%add legend & title
title('Accuracy for classification of nouns vs. verbs')
xlabel(sprintf('time in ms (center of %d ms sliding window)',cfg{plotcfg{1}(1)}.twidth*1000))
ax = gca;
L = [];
L(1) = plot(nan,nan,'.','MarkerEdgeColor',cmap(1,:),'MarkerSize',40,'Parent',ax);
L(2) = plot(nan,nan,'.','MarkerEdgeColor',cmap(2,:),'MarkerSize',40,'Parent',ax);
L(3) = plot(nan,nan,'.','MarkerEdgeColor',cmap(end,:),'MarkerSize',40,'Parent',ax);
L(4) = plot(nan,nan,'-k','LineWidth',5,'Parent',ax);
hl(end+1) = legend(L,{'accuracy noun vs verb','word position control','permuted class labels','significance',''},'Location','best');

%format
ylim([0.4 0.8])
xlim([-200 700])
for h = 1:length(hl)
    set(hl(h),'Linewidth',5)
end
set(hl(end),'FontSize',25)
set(ax,'FontSize',25)

print(gcf,fullfile(savedir,'plot3.eps'),'-depsc','-painters')

%% Plot X - different classifier methods (bayes vs svm vs blogreg)
plotcfg = {[6,16,17]};
taxis = round(time{plotcfg{1}(1)}*1000);

%%Stats
dat.dimord = 'chan_time';
dat.label = {'all'};
dat.time = time{plotcfg{1}(1)};
for nsub = 1:10
    acc_NV{nsub} = dat;
    acc_NV{nsub}.avg = acc{plotcfg{1}(1)}(:,nsub)';
    acc_pos{nsub} = dat;
    acc_pos{nsub}.avg = acc{plotcfg{1}(3)}(:,nsub)';  
end

cfgs = [];
cfgs.method = 'montecarlo';
cfgs.neighbours = [];
cfgs.statistic = 'ft_statfun_depsamplesT';
cfgs.correctm = 'cluster';
cfgs.correcttail = 'prob';
cfgs.clusteralpha = 0.05;
cfgs.clusterstatistic = 'maxsum';
cfgs.alpha = 0.025;
cfgs.numrandomization = 'all';
cfgs.design(1,1:2*nsub) = [ones(1,nsub) 2*ones(1,nsub)];
cfgs.design(2,1:2*nsub) = [1:nsub 1:nsub];
cfgs.ivar = 1;
cfgs.uvar = 2;
stat = ft_timelockstatistics(cfgs,acc_NV{:},acc_pos{:});

%save stats
save(fullfile(rootdir,'groupresults','sensor_stats_blogreg.mat'),'results','params')

%plot permutation distribution
maccshuf = mean(reshape(accshuf{plotcfg{1}(1)},length(taxis),[]),2);
accshufbounds = std(reshape(accshuf{plotcfg{1}(1)},length(taxis),[]),[],2);
figure('units','normalized','outerposition',[0 0 1 1]); hold on;
hl = [];
hl(1)= boundedline(taxis, maccshuf, accshufbounds, ...
    'alpha','cmap',cmap(end,:));

%plot data
for l = 1:length(plotcfg{1})
    datid = plotcfg{1}(l);
    
    macc = mean(acc{datid},2);
    accbounds = std(acc{datid},[],2);
    
    hl(1+l) = boundedline(taxis, macc, accbounds, ...
        'alpha','cmap',cmap(l,:));
    
end
%Plot significance indicator
idx_sig = double(stat.mask);
idx_sig(idx_sig==0) = NaN;
plot(taxis,(0.5-0.06)*idx_sig,'-k','LineWidth',5)

%add legend & title
title('Accuracy for classification of nouns vs. verbs')
xlabel(sprintf('time in ms (center of %d ms sliding window)',cfg{plotcfg{1}(1)}.twidth*1000))
ax = gca;
L = [];
L(1) = plot(nan,nan,'.','MarkerEdgeColor',cmap(1,:),'MarkerSize',40,'Parent',ax);
L(2) = plot(nan,nan,'.','MarkerEdgeColor',cmap(2,:),'MarkerSize',40,'Parent',ax);
L(3) = plot(nan,nan,'.','MarkerEdgeColor',cmap(3,:),'MarkerSize',40,'Parent',ax);
L(4) = plot(nan,nan,'.','MarkerEdgeColor',cmap(end,:),'MarkerSize',40,'Parent',ax);
L(5) = plot(nan,nan,'-k','LineWidth',5,'Parent',ax);

hl(end+1) = legend(L,{'Gaussian Naive Bayes','Support Vector Machines', 'Logistic Regression','permuted class labels','significance Logreg vs. GNB'},'Location','best');

%format
ylim([0.4 0.8])
xlim([-200 700])
for h = 1:length(hl)
    set(hl(h),'Linewidth',5)
end
set(hl(end),'FontSize',25)
set(ax,'FontSize',25)

print(gcf,fullfile(savedir,'plotSVMLogreg.eps'),'-depsc','-painters')

%% Plot4: decoding accuracy for NA/VA
plotcfg = {[14,15]};
taxis = round(time{plotcfg{1}(1)}*1000);

%plot permutation distribution
maccshuf = mean(reshape(accshuf{plotcfg{1}(1)},length(taxis),[]),2);
accshufbounds = std(reshape(accshuf{plotcfg{1}(1)},length(taxis),[]),[],2);
figure('units','normalized','outerposition',[0 0 1 1]); hold on;
hl = [];
hl(1)= boundedline(taxis, maccshuf, accshufbounds, ...
    'alpha','cmap',cmap(end,:));

%plot data
for l = 1:length(plotcfg{1})
    datid = plotcfg{1}(l);
    
    macc = mean(acc{datid},2);
    accbounds = std(acc{datid},[],2);
    
    hl(1+l) = boundedline(taxis, macc, accbounds, ...
        'alpha','cmap',cmap(l,:));
end

%add legend & title
title('Accuracy for classification of nouns vs. verbs')
xlabel(sprintf('time in ms (center of %d ms sliding window)',cfg{plotcfg{1}(l)}.twidth*1000))
ax = gca;
L = [];
L(1) = plot(nan,nan,'.','MarkerEdgeColor',cmap(1,:),'MarkerSize',40,'Parent',ax);
L(2) = plot(nan,nan,'.','MarkerEdgeColor',cmap(2,:),'MarkerSize',40,'Parent',ax);
L(3) = plot(nan,nan,'.','MarkerEdgeColor',cmap(end,:),'MarkerSize',40,'Parent',ax);
hl(end+1) = legend(L,{'Noun attached vs Verb attached','attachment label according to post-test','permuted class labels'},'Location','best');

%format
ylim([0.4 0.65])
xlim([-200 1900])
for h = 1:length(hl)
    set(hl(h),'Linewidth',5)
end
set(hl(end),'FontSize',25)
set(ax,'FontSize',25)

print(gcf,fullfile(savedir,'plot4_NAVA.eps'),'-depsc','-painters')

%% Generalize Nouns vs Verbs
clear all
savedir = fullfile('/project','3011210.01','MEG','figures','final');
rootdir = fullfile('/project','3011210.01','MEG','Classification');
load(fullfile(rootdir,'groupresults','NNVVFIN_generalize.mat'))

taxis1 = round(cfgcv.trainwind*1000);
taxis2 = round(cfgcv.time*1000);
[t1,t2,nrep,nsub] = size(accshuf);

%stats
a = cat(3,reshape(acc,[t1*t2,nsub]),permute(reshape(accshuf,[t1*t2,nrep,nsub]),[1,3,2]));
[results, params] = prevalenceCore(a);

macc = flipud(mean(acc,3));
figure;imagesc(macc)
caxis([0.4 0.6])
colorbar
set(gca, 'XTick',[1:round(t2/5):t2], 'XTickLabel', taxis2(1:round(t2/5):end))
set(gca, 'YTick',[1:round(t1/5):t1], 'YTickLabel', flip(taxis1(1:round(t1/5):end)))
title('Generalizing Noun-Verb Classifier to attachment labels.')
xlabel('Training time from onset final noun (in ms)')
ylabel('Testing time from onset Noun/Verb (in ms)')
print(gcf,fullfile(savedir,'plot_general.eps'),'-depsc','-painters')

%% Plot 5 - Searchlight approach in source space
clear all
rootdir = fullfile('/project','3011210.01','MEG','Classification');
savedir = fullfile('/project','3011210.01','MEG','figures','final');

cmap = brewermap(6,'Pastel1'); %for plotting lines
cmap_sig = brewermap(4,'OrRd'); %for plotting significance

fnames = '%s_acc_allparcels_new.mat';
load atlas_subparc374_8k
pindx = 1:length(atlas.parcellationlabel);
pindx([1 2 188 189]) = []; %ignore medial wall parcels

for nsub = 1:10
    subj = sprintf('sub-%.3d',nsub);
    load(sprintf(fullfile(rootdir,'%s','lcmv','searchlight',fnames),subj,subj)) 
    %map to full source space
    for p = 1:size(acc,2)  
        indx = pindx(p);
        a(:,indx,nsub) = acc(:,p);
        ashuf(:,:,indx,nsub) = accshuf(:,:,p);
    end
end
clear acc accshuf
a(a==0) = NaN;
ashuf(ashuf==0) = NaN;

%compute stats
[nt,nreps,nparc,nsubj] = size(ashuf);
atmp = reshape(cat(4,a,permute(ashuf,[1,3,4,2])),[nt*nparc,nsubj,nreps+1]);
[results, params] = prevalenceCore(atmp);
results = structfun(@(x)reshape(x,[nt nparc]),results,'UniformOutput',0);

%save stats
save(fullfile(rootdir,'groupresults','lcmv_searchlight_prevalence.mat'),'results','params')

%visualize cortical mesh
load(fullfile('/project/3011210.01/anatomy',subj,strcat(subj,'_sourcemodel.mat')));
source                = [];
source.brainordinate  = atlas;
source.brainordinate.pos = sourcemodel.pos;
source.label          = atlas.parcellationlabel;
source.time           = cfgcv.time;
source.dimord         = 'chan_time';
source.pow            = (mean(a,3).* (results.pcMN<0.05))';
source.mask           = results.pcMN<0.05;

%for exploration
cfgp                  = [];
cfgp.funparameter     = 'pow';
%cfgp.maskparameter    = source.mask;
ft_sourcemovie(cfgp, source);

%FIXME: save plot
%1. 60-160ms R18/R17
toi{1} = [nearest(cfgcv.time,0.060):nearest(cfgcv.time,0.160)];
roi{1} = find(strncmp(atlas.parcellationlabel,'R_18',4)); 
roi{1} = roi{1}(ismember(roi{1},find(all(results.pcMN(toi{1},:)<0.05))));
roiname{1} = 18;
viewplot{1} = [0.8 -10];
%2. 160-260 L37
toi{2} = [nearest(cfgcv.time,0.160):nearest(cfgcv.time,0.260)];
roi{2} = find(strncmp(atlas.parcellationlabel,'L_37',4));
roi{2} = roi{2}(ismember(roi{2},find(all(results.pcMN(toi{2},:)<0.05))));
roiname{2} = 37;
viewplot{2} = [-90 0];
%3. 340-440 L43
toi{3} = [nearest(cfgcv.time,0.340):nearest(cfgcv.time,0.440)];
roi{3} = find(strncmp(atlas.parcellationlabel,'L_45',4));
roi{3} = roi{3}(ismember(roi{3},find(all(results.pcMN(toi{3},:)<0.05))));
roiname{3} = 43;
viewplot{3} = [-90 0];

% create the 'upsampling matrix', to map parcels onto the vertices of the
% cortical sheet
x = zeros(0,1);
y = zeros(0,1);
for k = 1:numel(atlas.parcellationlabel)
    sel = atlas.parcellation==k;
    x = cat(1,x(:),find(sel));
    y = cat(1,y(:),ones(sum(sel),1)*k);
end
P = sparse(x,y,ones(numel(x),1),size(atlas.pos,1),numel(atlas.parcellationlabel));

%prepare for plotting
taxis = round(cfgcv.time*1000);
accmap = brewermap(20,'YlOrRd');
figure;
rgbplot(accmap)
hold on
colormap(accmap)
ax = gca;
ax.CLim = [0.50 0.70];
colorbar('Ticks',[0.5 0.55 0.6 0.65 0.7])
export_fig(fullfile(savedir,'colorbar_acc'),'-eps');

for i = 1:length(toi)
    shuf = reshape(permute(ashuf,[1,3,2,4]),length(taxis),nparc,[]);
    shuf = squeeze(nanmean(shuf(:,roi{i},:),2));
    dat = squeeze(nanmean(a(:,roi{i},:),2));
    %plot space
    mask = zeros(1,374);
    mask(roi{i}) = 1;
    figure('units','normalized','outerposition',[0 0 1 1]);
    ft_plot_mesh(atlas,'vertexcolor',P*mean(source.pow(:,toi{i}),2), ...
        'clim', [0.50 0.70],...
        'facealpha', P*double(any(results.pcMN(toi{i},:) < 0.05))', ...
        'colormap',brewermap(20,'YlOrRd'), ...
        'maskstyle', 'colormix', ...
        'contour',P*mask',...
        'contourcolor',[1 1 1]);
    lighting gouraud; material dull;
    view(viewplot{i});
    camlight;
    set(gcf,'color','w')
    title(sprintf('Average decoding accuracy from %d to %d ms for BA %d',taxis(toi{i}(1)),taxis(toi{i}(end)),roiname{i}));
    export_fig(fullfile(savedir,sprintf('lcmv_searchlight_NounVerb_BA%d',roiname{i})),'-png','-transparent','-m5');
    
    %plot time course
    maccshuf = mean(shuf,2);
    accshufbounds = std(shuf,[],2);
    figure('units','normalized','outerposition',[0 0 1 1]); hold on;
    hl = [];
    hl(1)= boundedline(taxis, maccshuf, accshufbounds, ...
        'alpha','cmap',cmap(end,:));
    
    macc = mean(dat,2);
    accbounds = std(dat,[],2);
    
    hl(2) = boundedline(taxis, macc, accbounds, ...
        'alpha','cmap',cmap(1,:));
    
    idx_sig = double(any(results.pcMN(:,roi{i})<0.05,2));
    idx_sig(idx_sig==0) = NaN;
    plot(taxis,(0.5-0.03)*idx_sig,'*','MarkerEdgeColor','k','MarkerSize',10)
    
    title(sprintf('Decoding accuracy Verb vs Noun at BA %i',roiname{i}))
    xlabel(sprintf('time in ms (center of %d ms sliding window)',cfgcv.twidth*1000))
    ax = gca;
    L = [];
    L(1) = plot(nan,nan,'.','MarkerEdgeColor',cmap(1,:),'MarkerSize',40,'Parent',ax);
    L(2) = plot(nan,nan,'.','MarkerEdgeColor',cmap(end,:),'MarkerSize',40,'Parent',ax);
    L(3) = plot(nan,nan,'*','MarkerEdgeColor','k','MarkerSize',20,'Parent',ax);
    hl(end+1) = legend(L,{'BA 18','permuted class labels','significance'},'Location','best');
    
    %format
    ylim([0.45 0.75])
    xlim([0 700])
    for h = 1:length(hl)
        set(hl(h),'Linewidth',5)
    end
    set(hl(end),'FontSize',25)
    set(ax,'FontSize',25)
    
    print(gcf,fullfile(savedir,sprintf('acc_lcmv_searchlight_parcel%i.eps',roiname{i})),'-depsc','-painters')
    
end

%% Plot - most informative sensors

clear all
rootdir = fullfile('/project','3011210.01','MEG','Classification');
savedir = fullfile('/project','3011210.01','MEG','figures','final');

cmap = brewermap(6,'Pastel1'); %for plotting lines
cmap_sig = brewermap(4,'OrRd'); %for plotting significance

fnames = {'nbayes_%s_lp01_20folds_1500feats_NNVVFIN_pca.mat'};

for nsub = 1:10
    subj = sprintf('sub-%.3d',nsub);
    filename = dir(sprintf(fullfile(rootdir,'%s','sensor',fnames{1}),subj,subj));
    load(fullfile(filename.folder,filename.name))
    %select output of interest (difference in normalised mean
    %between conditions)
    param = cell2mat(cfgcv.param(1:length(cfgcv.time)-1,1:cfgcv.nfolds));
    
    [nt, nfolds] = size(param);
    nfeat = size(param(1,1).Mudiff,2);
    
    Mudiff = reshape([param.Mudiff],nfeat,nt,nfolds);
    Mudiff = squeeze(mean(mean(reshape(Mudiff,60,31,nt,nfolds),4),2));
        
    Mudiff_chan{nsub} = Mudiff' * cfgcv.pca_unmixing;
end

for nsub = 1:10
    subj = sprintf('sub-%.3d',nsub);
    load(sprintf('/project/3011210.01/MEG/%s_dataclean',subj),'data')
    tmp.time  = cfgcv.time(1:end-1);
    tmp.label = data.label;
    tmp.avg   = Mudiff_chan{nsub}';
    tmp.dimord = 'chan_time';
    tmp.grad = data.grad;
    
    cfg = [];
    cfg.feedback = 'yes';
    cfg.method = 'template';
    cfg.neighbours = ft_prepare_neighbours(cfg,tmp);
    
    cfg.planarmethod = 'sincos';
    planar = ft_megplanar(cfg,tmp);
    
    cfg = [];
    planar = ft_combineplanar(cfg,planar);
    
    dat_planar{nsub} = planar;
end
clear data tmp

GA = ft_timelockgrandaverage([],dat_planar{:});

figure;
cfg = [];
cfg.layout = 'CTF275_helmet.mat';
ft_topoplotER(cfg,GA)


