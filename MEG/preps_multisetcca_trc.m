function [trc, tlck] = preps_multisetcca_trc(data, varargin)

output            = ft_getopt(varargin, 'output', 'rho');
dosmooth          = ft_getopt(varargin, 'dosmooth', 0);
output2           = ft_getopt(varargin, 'output2', 'average_mod');

switch output
    case 'rho'
        outputflag = 0;
    case 'Z'
        outputflag = 1;
    case 'Z_scaled'
        outputflag = 2;
end


if iscell(data)
    data = ft_appenddata([], data{:});
    
    p = cell(numel(data.label),1);
    for m = 1:numel(data.label)
        % assume evertything before the _ to denote a unique parcel
        tok = tokenize(data.label{m},'_');
        p{m} = tok{1};
    end
end

tlck = data;

if outputflag>0
    % temporarily replace 0 with their original nans, in order to properly
    % compute a d.f.
    for n = 1:length(tlck.trial)
        tlck.trial{n}(tlck.trial{n}==0) = nan;
    end
end

if dosmooth>0
    % do a boxcar smoothing of the time series
    for m = 1:size(tlck.trial,1)
        tlck.trial{m} = ft_preproc_smooth(squeeze(tlck.trial{m}),dosmooth); % use a smoothing kernel of odd number of samples
    end
end

nsmp = length(tlck.trial);
[nchan,nt] = size(tlck.trial{1});
for n = 1:length(tlck.trial)
    tlck.trial{n} = tlck.trial{n}(:,1:nt);
end

dat = permute(reshape(cell2mat(tlck.trial),nchan,nt,nsmp),[1,3,2]);
% subtract the mean across trials
dat = dat-nanmean(dat,2);

if outputflag>0
    % convert back to 0
    dof = squeeze(sum(isfinite(dat),2));
end
dat(~isfinite(dat)) = 0;

c = nan+zeros(size(dat,1),size(dat,1),size(dat,3));
for k = 1:numel(tlck.time{1})
    datx=dat(:,:,k);
    datc=datx*datx';
    c(:,:,k) = datc./sqrt(diag(datc)*diag(datc)');
    if outputflag>0
        % Fisher Z transform with standardization
        n = min(dof(:,k),dof(:,k)');
        c(:,:,k) = atanh(c(:,:,k))./sqrt(1./(n-3));
        c(~isfinite(c)) = 1; % replace the fisher z transformed infinity values with 1
    end
end

switch output2
    case 'average_mod'  
        % correction term assumes identity
        if outputflag<2
            trc.rho(:,:,:) = squeeze(mean(mean(c)))-1./nchan;
        else
            dat = c+repmat(diag(nan(nchan,1)),[1 1 size(dat,3)]);
            dat = reshape(dat,[],size(dat,3));
            %dat(~isfinite(dat)) = [];
            trc.rho(:,:,:) = nanmean(dat,1)./nanstd(dat,[],1);
        end
    case 'single_cross'
        trc.rho = reshape(c,[],numel(tlck.time{1}));
        for k = 1:nchan
            for m = 1:nchan
                label{k,m} = sprintf('%s_%s',tlck.label{k},tlck.label{m});
            end
        end
        trc.label = label(:);
end

trc.dimord   = 'chan_time';
trc.time     = tlck.time;