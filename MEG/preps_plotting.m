%% This script plots classification accuracy and saves figure as png

%% parameters
%which plot to generate
if ~exist('do_plotacc',     'var'), do_plotacc      = false;            end
if ~exist('do_plotgeneral', 'var'), do_plotgeneral  = false;            end
if ~exist('do_plotbl',      'var'), do_plotbl       = true;             end
if ~exist('do_confusion',   'var'), do_confusion    = false;            end
if ~exist('compute_freq',   'var'), compute_freq    = false;            end

if ~exist('file',       'var'), file        = '';               end

root_dir                = '/project/3011210.01/MEG';
[filepath, name, ext]   = fileparts(file);
parts                   = strsplit(name,'_');
mode                    = parts{1};
subj                    = parts{2};
save_file               = fullfile(root_dir,'figures',mode,strcat(name,'.png'));

%load in data
load(fullfile(root_dir,'Classification',subj,filepath, name))
%postprocess results
    if strcmp(cfgcv.mva,'ridgeregression_sa')
        [m,n] = size(stat);
        z     = size(stat{1},2);
        acc = reshape(vertcat(stat{:}),[m n z]);
        acc = cell2mat(squeeze(acc(:,:,1)));
        
        accshuf = reshape(vertcat(statshuf{:}),[m n z]);
        accshuf = cell2mat(squeeze(accshuf(:,:,1)));
    else
        c =  length(cfgcv.vocab);
        [m, n, ~] = size(stat);
        
        stat = cell2mat(stat);
        stat = struct2cell(stat);
        statshuf = cell2mat(statshuf);
        statshuf = struct2cell(statshuf);
        
        confusion = squeeze(stat(2,:,:));
        confusionshuf = squeeze(statshuf(2,:,:));
        mconf = mean(reshape(cat(3,confusion{:}),[c c m n]),4);
        stdconf = std(reshape(cat(3,confusion{:}),[c c m n]),[],4);
        
        acc = cell2mat(squeeze(stat(1,:,:)));
        accshuf = cell2mat(squeeze(statshuf(1,:,:)));
    end
    
if plot_featureimportance
    if strcmp(filepath,'lcmv')
        load atlas_subparc374_8k
        if strcmp(cfgcv.mva,'preps_naivebayes') %plotting normalized mean differences bt conditions
            %select output of interest (difference in normalised mean
            %between conditions)
            [t rep nfolds] = size(cfgcv.param);
            nfeat = size(cfgcv.param{3,2,20}.Mudiff,2)
            param = cell2mat(cfgcv.param);
            paramshuf = cell2mat(cfgcv.paramshuf);
            %reshape into TimeslicexRepxFoldxFeaturexTime
            
            Mudiff = reshape([param.Mudiff],nfeat,t,rep,nfolds);            
            MudiffShuf = reshape([paramshuf.Mudiff],nfeat,t,rep,nfolds);
            
            Mudiff = reshape(Mudiff,370,31,3,2,20);
            MudiffShuf = reshape(MudiffShuf,370,31,3,2,20);
            
            %reshape to full voxel space instead of parcels
            [nparc,nt,t,rep,nfold] = size(Mudiff);
            pindx = 1:length(atlas.parcellationlabel);
            pindx([1 2 188 189]) = []; %ignore medial wall parcels
            for p = 1:size(Mudiff,1)
                
                indx = pindx(p);
                m(indx,:,:,:,:) = Mudiff(p,:,:,:,:);
                mshuf(indx,:,:,:,:) = MudiffShuf(p,:,:,:,:);
            end
            clear Mudiff MudiffShuf
            %mean over reps,folds & time
            m = squeeze(mean(mean(mean(m,4),5),2));
            mshuf = squeeze(mean(mean(mean(mshuf,4),5),2));
            %visualize cortical mesh
            load(fullfile('/project/3011210.01/anatomy',subj,strcat(subj,'_sourcemodel.mat')));
            source                = [];
            source.brainordinate  = atlas;
            source.brainordinate.pos = sourcemodel.pos;
            source.label          = atlas.parcellationlabel;
            source.time           = cfgcv.time;
            source.dimord         = 'chan_time';
            source.pow            = squeeze(m(:,1,:)-mshuf(:,1,:));  
           
            cfgp                  = [];
            cfgp.funparameter     = 'pow';
            ft_sourcemovie(cfgp, source);
            %source.pow = squeeze(m(:,1,1)-mshuf(:,1,1));
            %source = ft_checkdata(source,'datatype','source');
            %figure;
            %ft_plot_mesh(source,'edgecolor','none','vertexcolor',source.pow);lighting gouraud;material dull;

        end
    else
        %plot topos
    end
end
%
if do_plotbl
    
    % compute confidence intervals
    accbounds       = std(acc,[],2);
    accshufbounds   = std(accshuf,[],2);
    macc            = mean(acc,2);
    maccshuf        = mean(accshuf,2);
    
    
    figure('units','normalized','outerposition',[0 0 1 1])
    hold on;
    taxis = round(cfgcv.time(1:nearest(cfgcv.time,cfgcv.time(end)-cfgcv.twidth))*1000);
    [hl, hp] = boundedline(taxis, [macc maccshuf], reshape([accbounds accshufbounds],[length(accbounds) 1 2]), ...
        'alpha');
    set(hl,'Linewidt',3);
    
    xlabel(sprintf('time in ms (center of %d ms sliding window)',cfgcv.twidth*1000))
    ylabel('classification accuracy')
    ylim([0.2 1])
    title(sprintf('classifier: %s - classes:%s - %s',cfgcv.mva, horzcat(cfgcv.vocab{:}),subj),'interpreter','none')
    legend('binary classes','permuted classes')
    ax = gca();
    ax.XAxisLocation = 'origin';
    ax.YAxisLocation = 'origin';
    ax.TickDir = 'out';
    box off;
    set(gca,'Layer','Top')
    set(legend,'Location','best')
    set(gca,'FontSize',25)
    set(gca,'LineWidth',4)
    
    export_fig(save_file,'-png');
end

%% plot frequency spectrum of decoding accuracy
if computefreq
    subjects = strsplit(sprintf('sub-%.3d ', [1:10]));
    subjects = subjects(~cellfun(@isempty, subjects));
    Fs  = 303;
    N = 303;
    t = ((1:N)/Fs)-0.2;
    for s = 1:10
        subj = subjects{s};
        file = sprintf('nbayes_%s_20folds_60feats_NNVVFIN_allt_cleaned.mat',subj);
        [filepath, name, ext]   = fileparts(file);
        load(fullfile(root_dir,'Classification',subj,name))
        stat = cell2mat(stat);
        stat = struct2cell(stat);
        acc = cell2mat(squeeze(stat(1,:,:)));
        macc            = mean(acc,2);
        fftacc= fft(macc,N)/(N/2);
        fftacc = fftacc(1:N/2+1);
        freqs = (Fs/N)*(1:N/2+1);
        spctrm(s,:) = real(fftacc).^2 + imag(fftacc).^2;
        %figure
        %plot(freqs,spctrm(s,:));xlim([0 60]);ylim([0 0.002]);
    end
    allspctrm = mean(spctrm);
    figure;
    hl = plot(freqs,allspctrm);xlim([0 60]);ylim([0 0.001]);
    
end



% %% plot classification accuracy
% if do_plotacc
% load(name)
%
% time = cfg.timeinfo-0.05;%fixme:outdated
% time = round(time*1000);
%
% figure('units','normalized','outerposition',[0 0 1 1])
% hold on;
% distributionPlot(accshuf','color','r');
% distributionPlot(acc','color','b');
% xticklabels(time)
% xlabel('time in ms (center of 100ms time slice)')
% ylabel('classification accuracy')
% ylim([0.3 1])
% title(sprintf('classifier: %s - classes:%s - %s',classifier, horzcat(classes{:}),subj),'interpreter','none')
% classstr = sprintf('%s vs. %s',classes{1},classes{2});
% h = get(gca,'Children');
% legend([h(3) h(15)],{classstr,'permuted'}')
% set(legend,'Location','best')
% set(gca,'FontSize',25)
% set(gca,'LineWidth',4)
%
% xt = xticks;
% if length(xt) > 10
%   xticks(xt(1:2:end))
%   time = time(1:2:end);
%   xticklabels(time)
% end
%
% fname = sprintf('%s/Figures/classacc_%s%s_%dfolds_%dfeats_%s%s',save_dir,subj,datasuffix,folds,numfeat,horzcat(classes{:}),suffix);
% export_fig(fname,'-png');
% clf;
% end
%
% %% plot general classification accuracy
% if do_plotgeneral
%     if ~exist('testtrig',   'var'), testtrig     = horzcat(trigger{6:7});   end
%
%     testpos = pos(cellfun(@(x) any(ismember(x,testtrig)),trigger));
%     name = fullfile(save_dir, subj, sprintf('classgeneral_%s_%dfeats_%sto%s',subj,numfeat,horzcat(classes{:}),horzcat(testpos{:})));
%     load(name)
%     load(strcat(name,'_shuf'))
%
%     figure('units','normalized','outerposition',[0 0 1 1])
%     h1 = plot(cfgtmp.timeinfo-0.05,accshuf,'color',[0,0,0]+0.5,'linewidt',4);
%     hold on
%     h2 = plot(cfgtmp.timeinfo-0.05,acc,'color','g','linewidt',4);
%     xlabel('time in s (center of 100ms time slice)')
%     ylabel('classification accuracy')
%     title(sprintf('classifier: %s - classes %s generalized to %s - %s',classifier, horzcat(classes{:}),horzcat(testpos{:}),subj),'interpreter','none')
%     classstr = sprintf('%s vs. %s',classes{1},classes{2});
%     legend([h1(1) h2],{'permuted',classstr}')
%     set(gca,'FontSize',25)
%     set(gca,'LineWidth',4)
%
%     traintim = cfgtmp.trainwindow*1000;
%     fname = sprintf('%s/Figures/classgeneral_%s_%dfolds_%dfeats_trainwindow%dto%d_%sto%s',...
%         save_dir,subj,folds,numfeat,traintim(1),...
%         traintim(2),horzcat(classes{:}),horzcat(testpos{:}));
%     export_fig(fname,'-png');
%     clf;
% end




