function spike = net_comp(windows)

%% Parameters
do_plot = 1;
alpha = 0.05;
bf = 0; % Bonferroni correction for power change analysis?
paired = 'paired';

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

if exist(out_folder,'dir') == 0
    mkdir(out_folder);
end

%% Get power
% Returns a structure with power for different time windows and spike vs
% not a spike, as well as indices with significant power change.
sig_dev = get_sd(alpha,bf);

%% Get network permanova statistics
perm_stats = get_perm;

%% Get network metrics
metrics = get_metrics;

%% Reduce perm stats and power to the times without significant power change
[sd_red,perm_red,metrics_red] = reduce_to_ns(sig_dev,perm_stats,metrics,paired);

%% Compare reduced network change in spike and not-spike
perm_red = assess_net_change(perm_red,alpha);

%% Compare reduced power change in spike and not-spike
sd_red = assess_power_change(sd_red,alpha);

%% Compare reduced metrics in spike and not-spike
metrics_red = assess_metric_change(metrics_red,alpha);


if do_plot
%% Plot network change
plot_net_change(perm_red,windows,out_folder,paired,1);
if strcmp(paired,'unpaired')
    plot_net_change(perm_red,windows,out_folder,paired,2);
end

%% Plot power change
plot_power_change(sd_red,windows,out_folder,paired,1);
if strcmp(paired,'unpaired')
    plot_power_change(sd_red,windows,out_folder,paired,2);
end

%% Plot metric change
plot_metric_change(metric_red,windows,out_folder,paired,1,'ge');

end

%% Classifier
spike = spike_classifier(perm_red);


end