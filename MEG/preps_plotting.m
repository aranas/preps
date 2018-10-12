%% This script plots classification accuracy and saves figure as png

%% parameters
%which plot to generate
if ~exist('do_plotacc',     'var'), do_plotacc      = false;            end
if ~exist('do_plotgeneral', 'var'), do_plotgeneral  = false;            end
if ~exist('do_plotbl',      'var'), do_plotbl       = false;            end
if ~exist('do_confusion',   'var'), do_confusion    = false;            end

if ~exist('filename',       'var'), filename        = '';               end

root_dir    = '/project/3011210.01/MEG';
fileparts   = strsplit(filename,'_');
mode        = fileparts{1};
subj        = fileparts{2};
save_file   = fullfile(root_dir,'figures',mode,strrep(strjoin(fileparts(2:end),'_'),'mat','png'));

%% plot classification accuracy
if do_plotacc
load(filename)

time = cfg.timeinfo-0.05;%fixme:outdated
time = round(time*1000);

figure('units','normalized','outerposition',[0 0 1 1])
hold on;
distributionPlot(accshuf','color','r');
distributionPlot(acc','color','b');
xticklabels(time)
xlabel('time in ms (center of 100ms time slice)')
ylabel('classification accuracy')
ylim([0.3 1])
title(sprintf('classifier: %s - classes:%s - %s',classifier, horzcat(classes{:}),subj),'interpreter','none')
classstr = sprintf('%s vs. %s',classes{1},classes{2});
h = get(gca,'Children');
legend([h(3) h(15)],{classstr,'permuted'}')
set(legend,'Location','best')
set(gca,'FontSize',25)
set(gca,'LineWidth',4)

xt = xticks;
if length(xt) > 10
  xticks(xt(1:2:end))
  time = time(1:2:end);
  xticklabels(time)
end

fname = sprintf('%s/Figures/classacc_%s%s_%dfolds_%dfeats_%s%s',save_dir,subj,datasuffix,folds,numfeat,horzcat(classes{:}),suffix);
export_fig(fname,'-png');
clf;
end

%% plot general classification accuracy 
if do_plotgeneral
    if ~exist('testtrig',   'var'), testtrig     = horzcat(trigger{6:7});   end

    testpos = pos(cellfun(@(x) any(ismember(x,testtrig)),trigger));
    filename = fullfile(save_dir, subj, sprintf('classgeneral_%s_%dfeats_%sto%s',subj,numfeat,horzcat(classes{:}),horzcat(testpos{:})));
    load(filename)
    load(strcat(filename,'_shuf'))
    
    figure('units','normalized','outerposition',[0 0 1 1])
    h1 = plot(cfgtmp.timeinfo-0.05,accshuf,'color',[0,0,0]+0.5,'linewidt',4);
    hold on
    h2 = plot(cfgtmp.timeinfo-0.05,acc,'color','g','linewidt',4);
    xlabel('time in s (center of 100ms time slice)')
    ylabel('classification accuracy')
    title(sprintf('classifier: %s - classes %s generalized to %s - %s',classifier, horzcat(classes{:}),horzcat(testpos{:}),subj),'interpreter','none')
    classstr = sprintf('%s vs. %s',classes{1},classes{2});
    legend([h1(1) h2],{'permuted',classstr}')
    set(gca,'FontSize',25)
    set(gca,'LineWidth',4)

    traintim = cfgtmp.trainwindow*1000;
    fname = sprintf('%s/Figures/classgeneral_%s_%dfolds_%dfeats_trainwindow%dto%d_%sto%s',...
        save_dir,subj,folds,numfeat,traintim(1),...
        traintim(2),horzcat(classes{:}),horzcat(testpos{:}));
    export_fig(fname,'-png');
    clf;
end

if do_plotbl
    load(fullfile(root_dir,'Classification',subj,filename))
    if strcmp(cfgcv.mva,'ridgeregression_sa')
        [m,n] = size(stat);
        acc = reshape(vertcat(stat{:}),[m n 4]);
        acc = cell2mat(squeeze(acc(:,:,1)));
        
        accshuf = reshape(vertcat(statshuf{:}),[m n 4]);
        accshuf = cell2mat(squeeze(accshuf(:,:,1)));
    else
        stat = cell2mat(stat);
        stat = struct2cell(stat);
        statshuf = cell2mat(statshuf);
        statshuf = struct2cell(statshuf);
        
        confusion = squeeze(stat(2,:,:));
        confusionshuf = squeeze(statshuf(2,:,:));
        mconf = mean(reshape(cat(3,confusion{:}),[5 5 19 50]),4);
        stdconf = std(reshape(cat(3,confusion{:}),[5 5 19 50]),[],4);
        
        acc = cell2mat(squeeze(stat(1,:,:)));
        accshuf = cell2mat(squeeze(statshuf(1,:,:)));
    end
    % compute confidence intervals
    accbounds       = std(acc,[],2);
    accshufbounds   = std(accshuf,[],2);
    macc            = mean(acc,2);
    maccshuf        = mean(accshuf,2);
    
    figure('units','normalized','outerposition',[0 0 1 1])
    hold on;
    taxis = round(cfgcv.time(1:nearest(cfgcv.time,cfgcv.time(end)-cfgcv.twidth))*1000);
    [hl, hp] = boundedline(taxis, [macc maccshuf], [accbounds accshufbounds],'alpha');
    set(hl,'linewidth',3)
    
    xticklabels(cfg.cvtime)
    xlabel(sprintf('time in ms (center of %d ms sliding window)',cfgcv.twidth*1000))
    ylabel('classification accuracy')
    ylim([0 1])
    title(sprintf('classifier: %s - classes:%s - %s',classifier, horzcat(classes{:}),subj),'interpreter','none')
    legend('multiclass','permuted classes')
    set(legend,'Location','best')
    set(gca,'FontSize',25)
    set(gca,'LineWidth',4)
    
    export_fig(save_file,'-png');
    clf;
end

if do_confusion
    
end


