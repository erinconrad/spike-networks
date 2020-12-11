function spike_comp

%{
WHY ARE THERE NANS IN GE???


do a more careful comparison of sd
%}

%% Parameters
alpha = 0.05;
%rm_rise = 1; 
met = 'ge';
windows = [0.1];
which_pre_rise = 0; % 2 is default
%comp_points = 2;  %3 is default
% 0 = absolute, 1 = z score, 2 = relative change from first one, 3 = like z
% score but subtracting first one

if which_pre_rise == 0
    wpr = 'manual_before_rise';
elseif which_pre_rise == 2
    wpr = 'cons';
end

%% Get file locations, load spike times and pt structure
locations = spike_network_files;
main_folder = locations.main_folder;
results_folder = [main_folder,'results/'];
script_folder = locations.script_folder;
addpath(genpath(script_folder));

bct_folder = locations.BCT;
addpath(genpath(bct_folder));
out_folder = [results_folder,'plots/'];

if exist(out_folder,'dir') == 0
    mkdir(out_folder);
end


%% Find windows for each spike in which no early spike rise
pre_spike = find_pre_spike_windows(windows);

%% Get signal power deviation
sig_dev = get_sd(alpha,0,0);

% convert this to be similar to pre_spike
pre_spike = convert_sd(sig_dev,windows,pre_spike);

%% Get network metrics
[metrics,is_spike_soz] = get_specified_metrics(windows,pre_spike,met,1);


%% Remove bad spikes
[metrics,is_spike_soz,pre_spike] = remove_bad_spikes(metrics,is_spike_soz,pre_spike);

%% Remove the time windows with an early spike rise, get slopes, and do significance testing
%metrics_red = remove_early_rise(metrics,pre_spike,wpr,comp_points,rm_rise,alpha);
include_times = include_which_times(metrics,met,pre_spike,nan);
metrics = generate_summary_stats(metrics,met,include_times);


%% Significance testing across patients
%metrics_red = agg_pts_test(metrics_red);

%% Number positive
%count_pts_sig(metrics_red,met,1,which_freq,1,0.05/7)

%% Build a classifier to predict spike vs not spike
%class_spike(metrics_red,met)

%% Figs
%scatter_all_pts(metrics_red,windows,met,which_pre_rise)
%make_fig13(metrics_red,windows,met,which_pre_rise)
%show_earliest_rise(metrics_red,windows,met,which_pre_rise)
new_rise_fig(metrics,met)

%% Compare metric between spikes in soz and spikes outside soz
%soz_comparison(metrics_red,is_spike_soz,met,out_folder)

%% Plot the avg in time windows across patients
%agg_pts_tw(metrics_red,met,windows,method,which_pre_rise)

%% Plot slopes across patients
%agg_pts_plot(metrics_red,met,windows,method)

%% Plot the distribution of slopes (single pt)
%plot_slopes(metrics_red,met,which_pt,windows)


end