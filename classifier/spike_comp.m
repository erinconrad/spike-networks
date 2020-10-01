function spike_comp

%{
Plan:
1) Look at each spike individually
2) Remove times of some rise above baseline
   - figure out what to do for patients where I can't get involved chs (75,
   83)
3) Calculate networks in each time period before that
4) Do same for non-spikes
5) Get slopes for node strength and ge in the pre-spike time for both
spikes and not spikes
6) Now I will have a bunch of slopes, 50 for spikes and 50 for not spikes
7) Compare slopes between spikes and not spikes; do classifier
%}

%% Parameters
met = 'ge';
windows = [0.1];
method = 'ttestpabs';
which_pt = 1;

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

if exist(out_folder,'dir') == 0
    mkdir(out_folder);
end


%% Find windows for each spike in which no early spike rise
pre_spike = find_pre_spike_windows(windows);

%% Get signal power deviation
sig_dev = get_sd(0.05,0,0);

%% Get network metrics
metrics = get_all_metrics(windows);

%% Remove the time windows with an early spike rise, get slopes, and do significance testing
metrics_red = remove_early_rise(metrics,pre_spike);

%% Significance testing across patients
metrics_red = agg_pts_test(metrics_red);

%% Plot slopes across patients
agg_pts_plot(metrics_red,met,windows,method)

%% Plot the distribution of slopes (single pt)
plot_slopes(metrics_red,met,which_pt,windows)


end