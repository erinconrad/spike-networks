function comp_soz

%% Parameters
alpha = 0.05;
rm_rise = 1; 
met = 'ers';
windows = [0.1];
which_pre_rise = 2; % 2 is default
comp_points = 3;  %3 is default
% 0 = absolute, 1 = z score, 2 = relative change from first one, 3 = like z
% score but subtracting first one

if which_pre_rise == 0
    wpr = 'manual_before_rise';
elseif which_pre_rise == 1
    wpr = 'before_rise';
    error('I do not do this anymore')
else
    wpr = 'cons';
end

%% Get file locations, load spike times and pt structure
locations = spike_network_files;
main_folder = locations.main_folder;
results_folder = [main_folder,'results/'];
data_folder = [main_folder,'data/'];
eeg_folder = [main_folder,'results/eeg_data/'];
script_folder = locations.script_folder;
addpath(genpath(script_folder));

bct_folder = locations.BCT;
addpath(genpath(bct_folder));
out_folder = [results_folder,'plots/'];
sig_dev_folder = [results_folder,'signal_deviation/manual/'];
perm_folder = [results_folder,'perm_stats/'];
ers_folder = [results_folder,'ers/'];
ns_folder = [results_folder,'metrics/manual/'];
adj_folder = [results_folder,'adj_mat/manual/'];
spike_rise_folder = [results_folder,'spike_rise/'];
pre_spike_folder = [results_folder,'/pre_spike/'];

%% Find windows for each spike in which no early spike rise
pre_spike = find_pre_spike_windows(windows);

%% Get signal power deviation
sig_dev = get_sd(alpha,0,0);
% convert this to be similar to pre_spike
pre_spike = convert_sd(sig_dev,windows,pre_spike);
%save([pre_spike_folder,'pre_spike.mat'],'pre_spike');

%% Get network metrics
[metrics,is_spike_soz] = get_specified_metrics(windows,pre_spike,met,1);

%% Remove the time windows with an early spike rise, get slopes, and do significance testing
metrics_red = remove_early_rise(metrics,pre_spike,wpr,comp_points,rm_rise,alpha);

%% Compare metric between spikes in soz and spikes outside soz
soz_comparison(metrics_red,is_spike_soz,met,out_folder)

end