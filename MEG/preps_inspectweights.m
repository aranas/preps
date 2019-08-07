%This script inspects the class means of the output when preps_decoding is
%called with "save_full" option
rootdir = '/project/3011210.01/MEG/Classification/all/mscca/';
file = dir(fullfile(rootdir,'*_NAVA_*_inspect_mscca3.mat'));

load(fullfile(rootdir,file.name))
nparcels = 370;
ntime = 31;


test = cell2mat(cfgcv.param);
fnames = fieldnames(test);
test = struct2cell(test);
test = squeeze(test(3,:,:,:));
[m,n,p] = size(test);
q = size(test{1},2);
test = reshape(vertcat(test{:}),[m,n,p,q]);
test = squeeze(nanmean(nanmean(test,3),2));

test = reshape(test,[m,nparcels,ntime]);
maxv = max(test(:));
minv = min(test(:));
for i = 1:size(test,1)
    figure;
    imagesc(squeeze(test(i,:,:))',[minv maxv])
    title(sprintf('100ms timeslice starting at %0.1f s',cfgcv.time(i)))
end

load atlas_subparc374_8k
pindx = 1:length(atlas.parcellationlabel);
pindx([1 2 188 189]) = []; %ignore medial wall parcels

mudiff = zeros(m,size(atlas.parcellationlabel,1),ntime);
for p = 1:nparcels
    indx = pindx(p);
    mudiff(:,indx,:) = test(:,p,:);
end

load('/project/3011210.01/anatomy/sub-001/sub-001_sourcemodel.mat');
source                = [];
source.brainordinate  = atlas;
source.brainordinate.pos = sourcemodel.pos;
source.label          = atlas.parcellationlabel;
source.dimord         = 'chan_time';
source.time           = 1:31;
source.pow            = squeeze(mudiff(1,:,:));

source = ft_checkdata(source,'datatype','source');

%plot both hemispheres simultaneously
pos = sourcemodel.pos;
n = 7842;
pos(1:n,1) = pos(1:n,1)+150;
pos(1:n,1:2) = -pos(1:n,1:2);
source.brainordinate.pos = pos;

cfgp                  = [];
cfgp.funparameter     = 'pow';
ft_sourcemovie(cfgp, source);
