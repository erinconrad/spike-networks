clear

%% Parameters
met = 'abs_power';
time_window = 0.1;

%% Get data and do basic processes
% Get data
orig = get_metric_data(met,time_window);

% Find windows for each spike in which no early spike rise
pre_spike = find_pre_spike_windows(orig.times);


% Find bad spikes
bad = identify_bad_spikes;

% Remove bad spikes
clean = remove_bad(orig,bad);




%% Do analyses

% Plot the metric over time for spikes and non spikes

