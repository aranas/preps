%preps_w2v_RSA
%loads stimuli with fastText info
%computes distance matrix between all word embeddings for selected
%words
%computes distance matrix between evoked activity of all selected words

clear all

%parameters
do_source = true;
datasuffix  = '_lp01';
doaverage = false;
%seltrig contains the selection of trigger values for both brain data and
%any amount of models to be tested, such that seltrig{n} = brain data and
%seltrig{1:n-1} = triggers for model comparison.
all_seltrig     = {[113,123,213,223],[115,125,215,225]...%verb mode, noun model
   [119,129,219,229]};
toverlap = 0.8;
twidth = 0.1;
begtim = 0;
endtim = 2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% For each subject compute & save %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% correlation between model RDM and neural RDM %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
jobid = {};
for nsub = 1:10
    subj = sprintf('sub-%.3d',nsub);%map to full source space
    suffix = '_overlap80';
    %% Load/select data & embeddings
    
    root_dir    = '/project/3011210.01/MEG';
    channelfile     = fullfile(root_dir,sprintf('%s_dataclean%s.mat',subj, datasuffix));
    lcmvfile        = fullfile(root_dir,strcat(subj,'_preps_lcmv_parc.mat'));
    load(channelfile,'data')
    
    for i = 1:length(all_seltrig)
        
        [seltrig{i}, pos{i}] = preps_help_collecttrig(subj, all_seltrig{i});
        sel = ones(length(data.trialinfo),1);
        sel             = sel & ismember(data.trialinfo(:,1),seltrig{i});
        cfg             = [];
        cfg.trials      = sel;
        datasel{i}         = ft_selectdata(cfg,data);
        
    end
    clear data
    upos = uniqueStrCell(pos{end});
    suffix = [strcat(upos{:}) suffix];
    
    %convert to source data
    if do_source
        parcel_indx = [];
        num_comp = 1;
        load(lcmvfile);
        source_parc.filterlabel = filterlabel;
        for i = 1:length(datasel)
            datasel{i} = preps_sensor2parcel(datasel{i},source_parc,num_comp,parcel_indx);
        end
        clear source_parc source filterlabel
    end
    
    if ~all(datasel{1}.trialinfo(:,2) == datasel{2}.trialinfo(:,2))
        warning('sentence ids in the different trial selections do NOT correspond! ')
        return
    end
    
    [brain_perword, model, all_pos, all_wids, all_embd] = select_embeddings(datasel,pos,doaverage);
    
    %% Create similarity matrices for data and Correlate with models  
    %neural RDM per parcel (searchlight across time) & correlate
    if ~exist('begtim','var'), begtim = datasel{end}.time{end}(1); end
    if ~exist('endtim','var'), endtim = datasel{end}.time{end}(end); end
    endtim = endtim - twidth;
    time = linspace(begtim, endtim, round(abs(begtim-endtim) ./ ...
        (twidth - toverlap * twidth)) + 1);
    
    for t = 1:length(time)
        btim = nearest(datasel{end}.time{1},time(t));
        etim = nearest(datasel{end}.time{1},time(t)+twidth);
        data = permute(brain_perword(:,:,btim:etim),[1,3,2]);
        
        jobid{nsub,t} = qsubfeval('preps_rsa_correlate',data,'euclidean',model,'spearman',[],1,...
            'memreq',(1024^3)*5,'timreq',60*30,'batchid',sprintf('preps_rsa_t%i_%s_%s',t,subj,suffix));
        
    end
    allinfo{nsub}.pos = all_pos;
    allinfo{nsub}.wid = all_wids;
    allinfo{nsub}.model = model;
end
%Collect results
for nsub = 1:10
    subj = sprintf('sub-%.3d',nsub);%map to full source space
 
    rho = [];
    for t = 1:length(time)
        load([jobid{nsub,t} '_output.mat'])
        out = argout{1};            
        for i = 1:length(out)
            rho(i,t,:) = out(i).rho;
        end
    end
    modelinfo = allinfo{nsub};
    save(fullfile(root_dir,'RSA','w2v',sprintf('%s_%s_rho.mat',subj,suffix)),'rho','modelinfo','time')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% Collect data of all subjects %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% compute significance & plot results %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Collect results
clear all
root_dir    = '/project/3011210.01/MEG';
load atlas_subparc374_8k

mapped_rho = [];
for nsub = 1:10
    subj = sprintf('sub-%.3d',nsub);%map to full source space
    load(fullfile(root_dir,'RSA','w2v',sprintf('%s_ADJANNVVFIN_overlap80_rho.mat',subj)))
    
    pindx = 1:length(atlas.parcellationlabel);
    pindx([1 2 188 189]) = []; %ignore medial wall parcels
    for i = 1:size(rho,1) %for every model
        for p = 1:size(rho,3)
            indx = pindx(p);
            mapped_rho(i,nsub,indx,:) = squeeze(rho(i,:,p));
        end
    end
end
for i = 1:size(mapped_rho,1)
avg_rho{i} = squeeze(mean(squeeze(mapped_rho(1,:,:,:))));
std_rho{i} = squeeze(std(squeeze(mapped_rho(1,:,:,:))));
end

mapped_rho(:,:,[1 2 188 189],:) = nan;
%% Stats
savedir = fullfile('/project','3011210.01','MEG','figures','final');
Nsub = 10;
modeli = 1;
load(fullfile('/project/3011210.01/anatomy',subj,strcat(subj,'_sourcemodel.mat')));
source                = [];
source.brainordinate  = atlas;
source.brainordinate.pos = sourcemodel.pos;
source.label          = atlas.parcellationlabel;
source.time           = time;
%source.dimord         = 'chan_time';


for nsub = 1:Nsub
    all_sources{nsub} = source;
    all_sources{nsub}.pow = squeeze(mapped_rho(modeli,nsub,:,:));
    all_bsl{nsub} = source;
    all_bsl{nsub}.pow = zeros(size(all_sources{nsub}.pow));
end

cfg = [];
cfg.method = 'montecarlo';
cfg.statistic = 'ft_statfun_depsamplesT';
cfg.parameter = 'pow';
cfg.correctm = 'cluster';
cfg.inside = true(1,length(atlas.parcellationlabel));
cfg.inside([1 2 188 189]) = false;
cfg.neighbours = 'biosemi16_neighb.mat'; %dummy neighbourhood structure
cfg.connectivity     = parcellation2connmat(atlas);
cfg.alpha = 0.05;
cfg.tail = 1;
cfg.design(1,1:2*Nsub) = [ones(1,Nsub) 2*ones(1,Nsub)];
cfg.design(2,1:2*Nsub) = [1:Nsub 1:Nsub];
cfg.numrandomization = 1000;
cfg.ivar = 1;
cfg.uvar = 2;
stat = ft_timelockstatistics(cfg,all_sources{:},all_bsl{:});


%% Visualize
%for exploration
source.pow = avg_rho{modeli};
source.std = std_rho{modeli};
source.mask = double(stat.posclusterslabelmat==1);

% cfgp                  = [];
% cfgp.funparameter     = 'pow';
% cfgp.maskparameter    = 'mask';
% ft_sourcemovie(cfgp, source);

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

rsamap = brewermap(20,'YlOrRd');
figure;
rgbplot(rsamap)
hold on
colormap(rsamap)
ax = gca;
ax.CLim = [0 0.01];
colorbar('Ticks',[0 0.005 0.01])
export_fig(fullfile(savedir,'colorbar_rsa'),'-eps');

taxis = round(time*1000);
[a,b] = find(stat.mask);
sig_tim = unique(b);
for i = sig_tim(1):3:sig_tim(end)
    tmpstat = stat.mask(:,i);
    roi = unique(cellfun(@(str) str(1:4), atlas.parcellationlabel(tmpstat), 'UniformOutput',false));
    figure('units','normalized','outerposition',[0 0 1 1]);
    ft_plot_mesh(atlas,'vertexcolor',P*source.pow(:,i), ...
        'clim', [0 0.01],...
        'facealpha', P*tmpstat, ...
        'colormap',rsamap, ...
        'maskstyle', 'colormix');
    lighting gouraud; material dull;
    view([-90 0]);
    camlight;
    set(gcf,'color','w')
    title(sprintf('Areas %s encoding word semantics around %i',horzcat(roi{:}),taxis(i+1)),'Interpreter','none');
    export_fig(fullfile(savedir,sprintf('RSA_semantics_%i',taxis(i+1))),'-png','-transparent','-m5');  
end

