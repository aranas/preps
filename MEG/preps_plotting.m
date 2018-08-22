%% This script plots classification accuracy and saves figure as png

%% Default parameters
if ~exist('subj',        'var'), subj         = 'pilot-005';                end
if ~exist('root_dir',    'var'), root_dir     = '/project/3011210.01/MEG/'; end
if ~exist('save_dir',    'var'), save_dir     = '/project/3011210.01/MEG/Classification'; end
if ~exist('classes',     'var'), classes      = {'ART', 'NN'};              end
if ~exist('classifier',  'var'), classifier   = 'preps_naivebayes';         end
if ~exist('folds',       'var'), folds        = 20;                         end
if ~exist('numfeat',     'var'), numfeat      = 250;                        end

filename = fullfile(save_dir, subj, sprintf('classacc_%s_%dfolds_%dfeats_%s',subj,folds,numfeat,horzcat(classes{:})));

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
title(sprintf('classifier: %s - classes:%s - %s',classifier, horzcat(classes{:}),subj))
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

fname = sprintf('%s/%s/classacc_%s_%dfolds_%dfeats_%s',save_dir,subj,subj,folds,numfeat,horzcat(classes{:}));
export_fig(fname,'-png');
clf;