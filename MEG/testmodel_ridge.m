% Y = X * Beta + noise;
%Let Y be word embeddings and X be brain data
%% channel meg data
subjectname = 'V1001';
%load in single trial data
tlck        = mous_db_getdata(subjectname,'meg_erf_allwords_02-nextword');

%select only words at least 500ms long
nsmp        = cellfun('size', tlck.trial, 2);
tlck        = ft_selectdata(tlck, 'rpt', find(nsmp>=211));
%Select sentence condition only
trign       = [1 5 2 6];
sel         = ismember(tlck.trialinfo(:,2),trign);

cfg         = [];
cfg.trials  = sel;
cfg.channel = {'all' '-EEG057','-EEG058'};
%cfg.latency = [-0.1 0.5]; does not work use ft_redefinetrial instead
tlck_sent   = ft_selectdata(cfg,tlck);
cfg = [];
cfg.toilim    = [-0.1 0.5];
tlck_sent = ft_redefinetrial(cfg,tlck_sent);

nsmp        = cellfun('size', tlck_sent.trial, 2);
tlck_sent        = ft_selectdata(tlck_sent, 'rpt', find(nsmp>=181));

clear tlck

[m,n] = size(tlck_sent.trial{1});
data = zeros(length(tlck_sent.trial),m,n);
for i = 1:length(tlck_sent.trial)
    data(i,:,:) = tlck_sent.trial{i};
end

datasel = data(:,:,1:31);
[numobs,n,p] = size(datasel);
datasel = reshape(datasel,[numobs,(n*p)]);
%% simulate betas and compute Y
q

y = NaN(1000,10);
x = datasel(1:1000,1:8000)*(10^13);
noise = 0.01*randn(1000,10);
beta = [ones(1000,10)*-1;ones(2000,10);ones(2000,10)*-1;ones(1500,10);ones(1500,10)*-1];
beta = beta*0.0025;

[x_scored, Mu, Sigma] = zscore(x);
x = randn(1000,8000);
y = x*beta;% + noise;
%split data in training and test set
x_train = x(1:800,:);
y_train = y(1:800,:);
x_test  = x(801:end,:);
y_test  = y(801:end,:);



[y_hat,beta_hat,lambda_hat,lambdas] = ridgeregression_sa(x_train,x_test,y_train,numobs/leaveout,10);%umobs/leaveout