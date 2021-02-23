function pretty_methods

%% Parameters
pt_name = 'HUP116';%'HUP082';
s = 50;
ch = [];
surround_time = 2;
windows = [-2:0.1:0];

freq_bands = [0.5 30;... %sub-gamma
    30 100;... % ;low gamma
    100 256]; % high_gamma
n_f = size(freq_bands,1);

%% Get file locations, load spike times and pt structure
locations = spike_network_files;
main_folder = locations.main_folder;
results_folder = [main_folder,'results/'];
eeg_folder = [results_folder,'eeg_data/'];
script_folder = locations.script_folder;
addpath(genpath(script_folder));

bct_folder = locations.BCT;
addpath(genpath(bct_folder));
out_folder = [results_folder,'plots/'];

if exist(out_folder,'dir') == 0
    mkdir(out_folder);
end

%% Load eeg data and pre-process
fname = [eeg_folder,pt_name,'_eeg.mat'];
spike = load(fname);
spike = spike.spike;
fs = spike(s).fs;
full_time = spike(s).surround_time;


ch = spike(s).biggest_dev;

data = spike(s).data;

% pre process
data = pre_processing(data,1,0,1,fs);
eeg = data;

% Get the specific eeg data
data = data(:,ch);

%% Get index windows
peak = round(size(data,1)/2);
n_windows = length(windows);
index_windows = zeros(n_windows,2);


            
for i = 1:n_windows
    index_windows(i,1) = peak + round(windows(i)*fs);
    index_windows(i,2) = peak + round(windows(i)*fs) + round((windows(2)-windows(1))*fs);
end

%% Abs power
% Get baseline
baseline = median(data); 

% get deviation from baseline (for all time points) and square to
% get power
dev = (abs((data - baseline)).^2); 
dev_windows = zeros(n_windows,1);
for t = 1:size(index_windows,1)
    dev_windows(t) = mean(dev(round(index_windows(t,1)):round(index_windows(t,2))));
end

%% Get ERS
X = data-baseline;
ers_windows = zeros(n_windows,n_f);
for t = 1:n_windows
    Xtemp = X(max(1,round(index_windows(t,1))):...
        min(length(X),round(index_windows(t,2))));
    powers = get_power(Xtemp,fs,freq_bands);
    ers_windows(t,:) = powers;
end

%% Plot the example spikes
if 1
figure
set(gcf,'position',[440 620 1000 176])
plot(linspace(-full_time,full_time,length(data)),data,'k','linewidth',2);
hold on
xlim([-surround_time surround_time]);
set(gca,'fontsize',20)
end

%% Add windows
if 0
figure
set(gcf,'position',[440 620 1000 176])
plot(linspace(-full_time,full_time,length(data)),data,'k','linewidth',2);
hold on
for t = 1:length(windows)
    plot([windows(t) windows(t)],ylim,'k--')
end
xlim([-surround_time 0]);
set(gca,'fontsize',20)
end

%% Show abs power
if 0
figure
set(gcf,'position',[440 620 1000 176])
plot([-surround_time:0.1:0],dev_windows,'k','linewidth',2);
hold on
for t = 1:length(windows)
    plot([windows(t) windows(t)],ylim,'k--')
end
end

%% Show ers
if 0
figure
set(gcf,'position',[440 620 1000 176])
plot([-surround_time:0.1:0],ers_windows(:,1),'k','linewidth',2);
hold on
plot([-surround_time:0.1:0],ers_windows(:,2),'k--','linewidth',2);
plot([-surround_time:0.1:0],ers_windows(:,3),'k-.','linewidth',2);
for t = 1:length(windows)
    plot([windows(t) windows(t)],ylim,'k--')
end
end

%% Show example window
if 0
ex_window = 19;
indices = [index_windows(ex_window,1):index_windows(ex_window,2)];
ex_data = data(indices);
figure
plot(linspace(windows(ex_window),windows(ex_window)+0.1,length(ex_data)),ex_data,'k','linewidth',2);
end

%% Show example spectrogram
if 0
figure
spectrogram(data(index_windows(1,1):index_windows(end,2)),128,120,128,fs,'yaxis');
caxis([-30 30])
end



end