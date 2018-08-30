%% This script plots classification accuracy and saves figure as png

%% Default parameters
if ~exist('subj',           'var'), subj         = 'pilot-005';                               end
if ~exist('root_dir',       'var'), root_dir     = '/project/3011210.01/MEG/';                end
if ~exist('save_dir',       'var'), save_dir     = '/project/3011210.01/MEG/Classification';  end
if ~exist('classes',        'var'), classes      = {'VA', 'NA'};                             end
if ~exist('classifier',     'var'), classifier   = 'preps_naivebayes';                        end
if ~exist('folds',          'var'), folds        = 20;                                        end
if ~exist('numfeat',        'var'), numfeat      = 250;                                       end
if ~exist('do_plotacc',     'var'), do_plotacc   = false;                                     end
if ~exist('do_plotgeneral', 'var'), do_plotgeneral   = false;                                 end
if ~exist('suffix',         'var'), suffix       = '';                                        end
pos = {'ART','NN','VVFIN','ADJA','APPR','NA','VA','Fill'};
trigger = {[111,114,121,124,211,214,221,224], %Determiner
    [112,115,122,125,212,215,222,225],        %Nouns
    [113,123,213,223],                        %Verbs
    [118,128,218,228],                        %Adjectives
    [116,126,216,226],                        %Preposition
    [219,229],                                %last word Noun attached
    [119,129],                                %last word Verb attached
    [30:39]};                                 %all words in filler sentences
%% plot classification accuracy
if do_plotacc
filename = fullfile(save_dir, subj, sprintf('classacc_%s_%dfolds_%dfeats_%s%s',subj,folds,numfeat,horzcat(classes{:}),suffix));

load(filename)
load(strcat(filename,'_shuf'))

time = cfg.timeinfo-0.05;
time = round(time*1000);

figure('units','normalized','outerposition',[0 0 1 1])
hold on;
distributionPlot(accshuf','color','r');
distributionPlot(acc','color','b');
xticklabels(time)
xlabel('time in ms (center of 100ms time slice)')
ylabel('classification accuracy')
ylim([0.2 0.9])
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

fname = sprintf('%s/Figures/classacc_%s_%dfolds_%dfeats_%s%s',save_dir,subj,folds,numfeat,horzcat(classes{:}),suffix);
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


