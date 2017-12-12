cfg = [];
cfg.dataset = strcat('/project/3011210.01/raw/301121001sopara11_1200hz_20171207_01.ds');

%hdr   = ft_read_header(cfg.dataset);

cfg.continuous      = 'yes';
cfg.channel         = {'EEG057', 'EEG058','EEG059','EEG060', 'EEG061','EEG062', 'EEG064'};
data2                = ft_preprocessing(cfg);

plot(data2.time{1}(1:5000),data2.trial{1}(1,1:5000),'r')
hold on
plot(data2.time{1}(1:5000),data2.trial{1}(2,1:5000),'b')
plot(data2.time{1}(1:5000),data2.trial{1}(3,1:5000),'g')


%cut data off at beginning and end
data_cut        = data;
data_cut.time   = data.time{1}(8400:144000);
data_cut.trial  = data.trial{1}(:,8400:144000);

%compute change in value at time point x (x - (x+1))
%for both triggers and photodiode timings
trig_change = data_cut.trial(1,1:(end-1)) - data_cut.trial(1,2:end);
trig_change(trig_change>0) = 0;
trig_change = abs(trig_change);
trig_indices = find(trig_change > 1);

phot_change = data_cut.trial(2,1:(end-1)) - data_cut.trial(2,2:end);
%phot_change(phot_change<2) = 0;
phot_change(phot_change>0) = 0;
phot_change = abs(phot_change);
phot_change(phot_change<2) = 0;

phot_indices = find(phot_change > 2.07);
phot_indices([75 94]) = [];
%throw away first photodiode reaction because no following trigger
phot_indices([1 2 10 11 19 20 28 29 37 38 46 47 55 56 64 65 73 74 82 83 91 92 100 101]) = [];

%distance between trigger should be same as distance between photodiode
%measurements
trig_timings    = data_cut.time(trig_indices);
phot_timings    = data_cut.time(phot_indices);
trig_dist       = trig_timings(2:end) - trig_timings(1:end-1); 
phot_dist       = phot_timings(2:end) - phot_timings(1:end-1); 

%test whether trigger timings coincide with duration times in logfile
trig_change2 = data.trial{1}(1,1:(end-1)) - data.trial{1}(1,2:end);
trig_change2(trig_change2>0) = 0;
trig_change2 = abs(trig_change2);
trig_indices2 = find(trig_change2 > 0);
trig_timings2    = data.time{1}(trig_indices2);
% correct for starting point and compare to logfile: (trig_timings2 - trig_timings2(1)) + 2.325