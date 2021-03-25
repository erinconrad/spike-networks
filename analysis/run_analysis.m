clear

%% Parameters
met = 'abs_power';
time_window = 0.1;


if strcmp(met,'abs_power')
    alpha = 0.05;
else
    alpha = 0.05/3;
end

%% Get data and do basic processes
% Get data
orig = get_metric_data(met,time_window);

% Find windows for each spike in which no early spike rise
pre_spike = find_pre_spike_windows(orig.times);

% Get the spike time powers
spike_power = get_spike_time_powers(orig);

% Find bad spikes
bad = identify_bad_spikes;

% Remove bad spikes
[clean,clean_spike,clean_pre] = remove_bad(bad,orig,spike_power,pre_spike);




%% Do analyses
if 0
% Plot the metric over time for spikes and non spikes
show_metric_over_time(clean,clean_spike,clean_pre,1,2)

% Get the number of patients with a pre-spike rise and an aggregate
% statistic

count_sig_power_rise(clean,clean_spike,clean_pre,1,alpha)


metric_rise_by_ch(clean,clean_spike,clean_pre,1,1,alpha)

agg_metric_rise_by_ch(clean,clean_spike,clean_pre,1,alpha)


end

