%load all subjects
for sub = 1:10
load(sprintf('/project/3011210.01/MEG/Classification/sub-%.3d/sensor/regress_sub-%.3d_8folds_lambda1_ADJANNVVFIN_w2vdim300_eval2_w2v.mat',sub,sub))
 
[m,n] = size(stat);
z = size(stat{1}{1},2);
x = vertcat(stat{:});
rho = reshape(x(:,1),[m n]);
rho = reshape(cell2mat(rho),[m z n]);
rhoall(sub,:,:) = mean(rho,3);
bounds(sub,:,:) = std(rho,[],3);

x = vertcat(statshuf{:});
rhoshuf = reshape(x(:,1),[m n]);
rhoshuf = reshape(cell2mat(rhoshuf),[m z n]);
rhoshufall(sub,:,:,:) = rhoshuf;
shufbounds(sub,:,:) = std(rhoshuf,[],3);
end

rho = squeeze(mean(rhoall));
rhoshuf = squeeze(mean(reshape(permute(rhoshufall,[1 4 2 3]),[500 44 300])));
bound = squeeze(std(rhoall));
boundshuf = squeeze(std(reshape(permute(rhoshufall,[1 4 2 3]),[500 44 300])));

%check which dimensions are most consistent across subjects
%[a, b] = sort(mean(bound),'ascend');
[a, b] = sort(mean(rho),'descend');
rho = rho(:,b(1:5));
rhoshuf = rhoshuf(:,b(1:5));
bound = bound(:,b(1:5));
boundshuf = boundshuf(:,b(1:5));


taxis = round(cfgcv.time(1:nearest(cfgcv.time,cfgcv.time(end)-cfgcv.twidth/2))*1000);
figure('units','normalized','outerposition',[0 0 1 1])
hold on;
[hl, hp] = boundedline(taxis, rhoshuf, reshape(boundshuf,[length(boundshuf) 1 5]),'k','alpha');
[hl2, hp2] = boundedline(taxis, rho, reshape(bound,[length(bound) 1 5]),'alpha');
set(hl2,'Linewidt',3);
 
xlabel(sprintf('time in ms (center of %d ms sliding window)',cfgcv.twidth*1000))
ylabel('rho of prediction with true value')
title('mean correlation of true & predicted fastText word embedding across subjects for 5 dimensions (L2 ridge)')
ax = gca();
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
ax.TickDir = 'out';
box off;
set(gca,'Layer','Top')
set(legend,'Location','best')
set(gca,'FontSize',25)
set(gca,'LineWidth',4)
    
export_fig('/project/3011210.01/MEG/figures/regress/regress_8folds_lambda1_ADJANNVVFINeval2_w2v.mat','-png');


%evaluate on all dimensions
%take 5 highest correlations - are those consistent across subjects?
load(sprintf('/project/3011210.01/MEG/Classification/sub-%.3d/sensor/regress_sub-%.3d_8folds_lambda1_ADJANNVVFIN_w2vdim300_eval2_w2v.mat',sub,sub))
[m,n] = size(stat);
z = size(stat{1}{1},2);
x = vertcat(stat{:});
rho = reshape(x(:,1),[m n]);
rho = reshape(cell2mat(rho),[m z n]);

x = vertcat(statshuf{:});
rhoshuf = reshape(x(:,1),[m n]);
rhoshuf = reshape(cell2mat(rhoshuf),[m z n]);

[y i] = sort(squeeze(mean(mean(rho(11:44,:,:)),3)),'descend');

figure;
hold on
plot(taxis,squeeze(mean(rhoshuf,3)),'k')
plot(taxis,squeeze(mean(rho(:,i(1:5),:),3)))