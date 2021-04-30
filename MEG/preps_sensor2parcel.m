function [parceldata] = preps_sensor2parcel(data, source,ncomp,varargin)

% this function projects the sensor-level data into source space
% ensure the channel labels in the data to match the order of the channels
% in the spatial filter, and compute the parcel specific time courses
[a,b]      = match_str(source.filterlabel, data.label);
parceldata = data;
parceldata.trial = cell(size(data.trial));
if isempty(varargin{1})
    for parcel_indx = 1:length(source.label)
        
        %indx       = 1:min(5,size(source.F{parcel_indx},1));
        parceldata.trial = cellfun(@(x,y) [x;y], parceldata.trial, source.F{parcel_indx}(1:ncomp,a)*cellrowselect(data.trial,b),'UniformOutput', false);
        
    end
    ncomp = parcel_indx;
    parceldata.label   = source.label(1:ncomp);
else
    parcel_indx = varargin{1};
    indx       = 1:min(ncomp,size(source.F{parcel_indx},1));
    parceldata.trial = source.F{parcel_indx}(indx,a)*cellrowselect(data.trial,b);
    ncomp = size(indx,2);
    parceldata.label   = source.label(parcel_indx);
end

cfg                = [];
cfg.demean         = 'yes';
cfg.baselinewindow = [-0.2 0];
parceldata         = ft_preprocessing(cfg, parceldata);
parceldata         = removefields(parceldata, {'grad' 'elec'});