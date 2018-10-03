function [trl, event] = ft_trialfun_preps(cfg)

% FT_TRIALFUN_GENERAL determines trials/segments in the data that are
% interesting for analysis, using the general event structure returned
% by read_event. This function is independent of the dataformat
%
% The trialdef structure can contain the following specifications
%   cfg.trialdef.eventtype  = string
%   cfg.trialdef.eventvalue = number, string or list with numbers or strings
%   cfg.trialdef.prestim    = latency in seconds (optional)
%   cfg.trialdef.poststim   = latency in seconds (optional)
%
% If you want to read all data from a continous file in segments, you can specify
%    cfg.trialdef.triallength = duration in seconds (can be Inf)
%    cfg.trialdef.ntrials     = number of trials
%
% If you specify
%   cfg.trialdef.eventtype  = '?'
% a list with the events in your datafile will be displayed on screen.
%
% If you specify
%   cfg.trialdef.eventtype = 'gui'
% a graphical user interface will allow you to select events of interest.
%
% See also FT_DEFINETRIAL, FT_PREPROCESSING

% Copyright (C) 2005-2018, Robert Oostenveld
%
% This file is part of FieldTrip, see http://www.fieldtriptoolbox.org
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id$

% some events do not require the specification a type, pre or poststim period
% in that case it is more convenient not to have them, instead of making them empty
if ~isfield(cfg, 'trialdef')
  cfg.trialdef = [];
end
if isfield(cfg.trialdef, 'eventvalue')  && isempty(cfg.trialdef.eventvalue   ), cfg.trialdef = rmfield(cfg.trialdef, 'eventvalue' ); end
if isfield(cfg.trialdef, 'prestim')     && isempty(cfg.trialdef.prestim      ), cfg.trialdef = rmfield(cfg.trialdef, 'prestim'    ); end
if isfield(cfg.trialdef, 'poststim')    && isempty(cfg.trialdef.poststim     ), cfg.trialdef = rmfield(cfg.trialdef, 'poststim'   ); end
if isfield(cfg.trialdef, 'ntrials')     && isempty(cfg.trialdef.ntrials      ), cfg.trialdef = rmfield(cfg.trialdef, 'ntrials'    ); end

% default file formats
cfg.eventformat   = ft_getopt(cfg, 'eventformat');
cfg.headerformat  = ft_getopt(cfg, 'headerformat');
cfg.dataformat    = ft_getopt(cfg, 'dataformat');

% get the header, among others for the sampling frequency
if isfield(cfg, 'hdr')
  ft_info('using the header from the configuration structure\n');
  hdr = cfg.hdr;
else
  % read the header, contains the sampling frequency
  ft_info('reading the header from ''%s''\n', cfg.headerfile);
  hdr = ft_read_header(cfg.headerfile, 'headerformat', cfg.headerformat);
end

% get the events
if isfield(cfg, 'event')
  ft_info('using the events from the configuration structure\n');
  event = cfg.event;
else
  ft_info('reading the events from ''%s''\n', cfg.headerfile);
  event = ft_read_event(cfg.headerfile, 'headerformat', cfg.headerformat, 'eventformat', cfg.eventformat, 'dataformat', cfg.dataformat);
end

trl = [];
val = [];

% start by selecting all events
sel = true(1, length(event)); % this should be a row vector

% select all events of the specified type
if isfield(cfg.trialdef, 'eventtype') && ~isempty(cfg.trialdef.eventtype)
  for i=1:numel(event)
    sel(i) = sel(i) && ismatch(event(i).type, cfg.trialdef.eventtype);
  end
elseif ~isfield(cfg.trialdef, 'eventtype') || isempty(cfg.trialdef.eventtype)
  % search for trial events
  for i=1:numel(event)
    sel(i) = sel(i) && ismatch(event(i).type, 'trial');
  end
end

% select all events with the specified value
if isfield(cfg.trialdef, 'eventvalue') && ~isempty(cfg.trialdef.eventvalue)
  for i=1:numel(event)
    sel(i) = sel(i) && ismatch(event(i).value, cfg.trialdef.eventvalue);
  end
end

% convert from boolean vector into a list of indices
sel = find(sel);

for i=sel
  % catch empty fields in the event table and interpret them meaningfully
  if isempty(event(i).offset)
    % time axis has no offset relative to the event
    event(i).offset = 0;
  end
  if isempty(event(i).duration)
    % the event does not specify a duration
    event(i).duration = 0;
  end
  % determine where the trial starts with respect to the event
  if ~isfield(cfg.trialdef, 'prestim')
    trloff = event(i).offset;
    trlbeg = event(i).sample;
  else
    % override the offset of the event
    trloff = round(-cfg.trialdef.prestim*hdr.Fs);
    % also shift the begin sample with the specified amount
    trlbeg = event(i).sample + trloff;
  end
  % determine the number of samples that has to be read (excluding the begin sample)
  if ~isfield(cfg.trialdef, 'poststim')
    trldur = max(event(i).duration - 1, 0);
  else
    % this will not work if prestim was not defined, the code will then crash
    trldur = round((cfg.trialdef.poststim+cfg.trialdef.prestim)*hdr.Fs) - 1;
  end
  trlend = trlbeg + trldur;
  % add the beginsample, endsample and offset of this trial to the list
  % if all samples are in the dataset
  if trlbeg>0 && trlend<=hdr.nSamples*hdr.nTrials
    trl = [trl; [trlbeg trlend trloff]];
    if isnumeric(event(i).value)
      val = [val; event(i).value];
    elseif ischar(event(i).value) && numel(event(i).value)>1 && (event(i).value(1)=='S'|| event(i).value(1)=='R')
      % on brainvision these are called 'S  1' for stimuli or 'R  1' for responses
      val = [val; str2double(event(i).value(2:end))];
    else
      val = [val; nan];
    end
  end
end

% append the vector with values
if ~isempty(val) && ~all(isnan(val)) && size(trl,1)==size(val,1)
  trl = [trl val];
end

%% add ID using logfile
load('/home/language/sopara/Prepositionalphrases/preps/Stimuli/preps_stimuli.mat')

fid                  = fopen(cfg.logfile);
format               = ['%f %f %s %f %f %f %s %f %f','%*[^\n]'];
logtxt               = textscan(fid,format,'Headerlines',6,'Delimiter','\t');
rowHeadings          = {'Subject','Trial','Word','Condition','PairNum'...
    'VerbNum','Attachment','Time','Duration'};
logtxt               = cell2struct(logtxt,rowHeadings,2);

%correct erroneous trigger from trl file
if logtxt.Subject(1)==3
    trl(2442,:) = [];
 elseif logtxt.Subject(1)==6
     trl(2351,:) = [];
elseif logtxt.Subject(1)==9
    trl(2089,:) = [];
elseif logtxt.Subject(1)==10
    trl([1033 2344],:) = [];
end


if size(logtxt.Word,1) ~= size(trl,1)
    warning('logfile does not correspond to trigger numbers!!!')
end

for i = 1:length(logtxt.Word)
logtxt.Word = strrep(logtxt.Word,'.','');
end

for i = 1:length(stimuli)
    ind2            = [];
    ind         = find(logtxt.Condition==stimuli(i).condition);
    ind2        = logtxt.PairNum(ind)==stimuli(i).pair_num;
    ind         = ind(ind2);
    if length(ind) > 9
        [~,index] = sort(trl(ind,4));
        if strcmp(stimuli(i).attachment,'NA');
            ind = ind(index(1:9));
        elseif strcmp(stimuli(i).attachment,'VA');
            ind = ind(index(end-8:end));
        end
    end
    if ~strcmp(stimuli(i).question,'')
        ind(end+1)        = ind(end)+1;
    end
    trl(ind,5)        = stimuli(i).id;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION returns true if x is a member of array y, regardless of the class of x and y
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function s = ismatch(x, y)
if isempty(x) || isempty(y)
  s = false;
elseif ischar(x) && ischar(y)
  s = strcmp(x, y);
elseif isnumeric(x) && isnumeric(y)
  s = ismember(x, y);
elseif ischar(x) && iscell(y)
  y = y(strcmp(class(x), cellfun(@class, y, 'UniformOutput', false)));
  s = ismember(x, y);
elseif isnumeric(x) && iscell(y) && all(cellfun(@isnumeric, y))
  s = false;
  for i=1:numel(y)
    s = s || ismember(x, y{i});
  end
else
  s = false;
end

