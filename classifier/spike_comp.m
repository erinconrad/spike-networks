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


do a more careful comparison of sd
%}

%% Parameters
rm_rise = 1;
met = 'F';
windows = [0.1];
method = 'ttestp';
which_pt = 1;
which_pre_rise = 2;
comp_points = 3;  % 0 = absolute, 1 = z score, 2 = relative change from first one

if which_pre_rise == 0
    wpr = 'manual_before_rise';
elseif which_pre_rise == 1
    wpr = 'before_rise';
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

if exist(out_folder,'dir') == 0
    mkdir(out_folder);
end


if exist(pre_spike_folder,'dir') == 0
    mkdir(pre_spike_folder);
end


%% Find windows for each spike in which no early spike rise
pre_spike = find_pre_spike_windows(windows);

%% Get signal power deviation
sig_dev = get_sd(0.05,0,0);
% convert this to be similar to pre_spike
pre_spike = convert_sd(sig_dev,windows,pre_spike);
save([pre_spike_folder,'pre_spike.mat'],'pre_spike');

%% Get network metrics
metrics = get_all_metrics(windows,pre_spike);

%% Remove the time windows with an early spike rise, get slopes, and do significance testing
metrics_red = remove_early_rise(metrics,pre_spike,wpr,comp_points,rm_rise);

%% Significance testing across patients
metrics_red = agg_pts_test(metrics_red);


%% Plot the avg in time windows across patients
agg_pts_tw(metrics_red,met,windows,method)

%% Plot slopes across patients
agg_pts_plot(metrics_red,met,windows,method)

%% Plot the distribution of slopes (single pt)
%plot_slopes(metrics_red,met,which_pt,windows)


end