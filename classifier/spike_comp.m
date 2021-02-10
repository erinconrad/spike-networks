%function spike_comp

%{
WHY ARE THERE NANS IN GE???
confirm pre spike power rise, seems quite high quite early
why ers not significant for 0.2 s windows
check everything
post-spike
think more about auc measure = why so sig

manually get biggest dev ch  =  need this because by automatically finding
it and then taking biggest for spike and average for not spike, I am
biasing myself toward finding bigger earlier changes for spike
%}

%% Clear
clear

%% Parameters
do_cumulative = 0;
do_plot = 0;
%rm_rise = 1; 
met = 'ers';
windows = [0.1];
rm_rise = 1;
%which_pre_rise = 0; % 2 is default
%comp_points = 2;  %3 is default
% 0 = absolute, 1 = z score, 2 = relative change from first one, 3 = like z
% score but subtracting first one

%{
if which_pre_rise == 0
    wpr = 'manual_before_rise';
elseif which_pre_rise == 2
    wpr = 'cons';
end
%}

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
pre_spike = multi_reviewer_pre_spike(windows);

%% Compare two reviewers
earliest_rise = compare_two_reviewers(pre_spike);

%% Get signal power deviation
sig_dev = get_sd(0.05,0,0);

% convert this to be similar to pre_spike
pre_spike = convert_sd(sig_dev,windows,pre_spike);

%% Get network metrics
if strcmp(met,'sd')
    [metrics,is_spike_soz] = get_specified_metrics(windows,pre_spike,'ers',1);
else
    [metrics,is_spike_soz] = get_specified_metrics(windows,pre_spike,met,1);
end


%% Remove bad spikes
[metrics,is_spike_soz,pre_spike] = remove_bad_spikes(metrics,is_spike_soz,pre_spike);

%% Remove the time windows with an early spike rise, get slopes, and do significance testing
include_times = include_which_times(metrics,met,pre_spike,nan);
metrics = generate_summary_stats(metrics,met,include_times,rm_rise,is_spike_soz,do_cumulative);

%% Build a classifier to predict spike vs not spike
%tbl = class_spike(metrics,met);

%% Figs
plot_auc(metrics,met,out_folder,do_plot)
plot_short(metrics,met,3,earliest_rise,out_folder,do_plot);
soz_comparison(metrics,met,out_folder)




%end