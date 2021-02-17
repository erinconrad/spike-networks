%function spike_comp

%{
post-spike

look at uninvolved chs

%}

%% Clear
clear

%% Parameters
do_auto = 1;
do_cumulative = 0;
do_plot = 0;
%rm_rise = 1; 
met = 'ns_avg';
windows = 0.1;
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

if do_auto
    if strcmp(met,'ns_big')
        met = 'ns_auto';
    elseif strcmp(met,'ns_avg')
        met = 'ns_avg';
    else
        met = [met,'_auto'];
    end
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
pre_spike = multi_reviewer_pre_spike(windows);

%% Compare two reviewers
[earliest_rise,pt_rise] = compare_two_reviewers(pre_spike);

%% Get signal power deviation
sig_dev = get_sd;

% convert this to be similar to pre_spike
pre_spike = convert_sd(sig_dev,windows,pre_spike,met);

%% Get network metrics
[metrics,is_spike_soz] = get_specified_metrics(windows,pre_spike,met);


%% Remove bad spikes
[metrics,is_spike_soz,pre_spike,pt_rise] = remove_bad_spikes(metrics,is_spike_soz,pre_spike,pt_rise);

%% Remove the time windows with an early spike rise, get slopes, and do significance testing
include_times = include_which_times(metrics,met,pre_spike,nan);
metrics = generate_summary_stats(metrics,met,include_times,rm_rise,is_spike_soz,do_cumulative);

%% Build a classifier to predict spike vs not spike
%tbl = class_spike(metrics,met);

%% Sanity checks
% sanity_checks(metrics,met,pt_rise)

%% Figs
plot_auc(metrics,met,out_folder,do_plot);
plot_short_both(metrics,met,1,earliest_rise,out_folder,do_plot,rm_rise);
%soz_comparison(metrics,met,out_folder)


%plot_short(metrics,met,1,earliest_rise,out_folder,do_plot);




%end