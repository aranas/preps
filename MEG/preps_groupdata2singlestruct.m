function out = preps_groupdata2singlestruct(groupdata, subjects)

T = zeros(numel(groupdata{1}.trial),numel(groupdata{1})+1);
for s = 1:numel(groupdata)
  label = cell(numel(groupdata{s}.label),1);
  for m = 1:numel(label)
    label{m,1} = sprintf('%s_chan%03d',subjects{s},m);
  end
  groupdata{s}.label=label;
  T(:,s) = groupdata{s}.trialinfo(:,end); %adding IDs as presented to each subject
end
T = [T groupdata{s}.trialinfo(:,1:2)]; %adding IDs as aligned
out = groupdata{1};
label = out.label;
for s = 2:numel(groupdata)
  out.trial = cellcat(1,out.trial,groupdata{s}.trial);
  label = cat(1,label,groupdata{s}.label);
end
out.label     = label;
