function [w2v,data,labels] = preps_getw2v(data,varargin)

ncluster             = ft_getopt(varargin, 'ncluster', 1); %if 0 no clustering is done, if >0 than x number of clusters are computed

load preps_stimuli
fprintf('retrieving continuous word embedding values\n')
indsel = false(1,size(data.trial,2));
for smp = 1:size(data.trial,2)
    id          = data.trialinfo(smp,2);
    trig        = num2str(data.trialinfo(smp,1));
    w2vtmp      = stimuli(id).words(str2num(trig(end))).w2v;
    if size(w2vtmp,2)~=0
        w2v.feat(smp,:)     = stimuli(id).words(str2num(trig(end))).w2v;
        w2v.pos{smp}        = stimuli(id).words(str2num(trig(end))).pos{1};
        w2v.word{smp}       = stimuli(id).words(str2num(trig(end))).word{1};
        w2v.id{smp}         = id;
        indsel(smp) = 1;
    end
end
%only keep trials for which w2v is known
w2v.feat        = w2v.feat(indsel,:);
w2v.pos         = w2v.pos(indsel);
w2v.word        = w2v.word(indsel);
labels          = w2v.feat;

cfg             = [];
cfg.trials      = indsel;
data            = ft_selectdata(cfg,data);

if ncluster>1
    fprintf('converting word embeddings into categorical labels through spatial clustering')
    [C,IA,IC]       = unique(w2v.feat,'rows'); % get rid of repeats
    embeddings      = w2v.feat(IA,:);
    
    [~, embeddings] = pca(embeddings',50);
    embeddings      = zscore(embeddings);
    
    %spatial clustering
    [cidx, cmeans] = kmeans(embeddings,ncluster,'dist','sqeuclidean');
    % [silh,h] = silhouette(embeddings,cidx,'sqeuclidean');
    labels = cidx(IC);
    
    %t-sne to visually inspect natural clustering
    % rng default
    % Y3 = tsne(embeddings,cidx,3,50);
    % figure
    % scatter3(Y3(:,1),Y3(:,2),Y3(:,3),30,cidx,'filled');
    % %view(-93,14)
    %
    % text(Y3(:,1),Y3(:,2),Y3(:,3),w2v.word,'Fontsize',10);
    %how many verb-noun pairs end up in the same cluster?
    data.trialinfo(:,3) = cidx(IC);
    for id = 1:200
        indx = find(data.trialinfo(:,2)==id);
        if length(indx) == 4
            same(id) = data.trialinfo(indx(2),3) == data.trialinfo(indx(3),3);
        else
            id
        end
    end
end