function include_times = include_which_times(metrics,met,pre_spike,start_sec)

%% Get the time windows of the metric
times = metrics.time.freq(1).(met).pt(1).times;

%% If start_sec is nan, reset to the first time
if isnan(start_sec)
    start_sec = times(1); 
end

% Loop through patients
for p = 1:length(pre_spike)
    
    %% Identify the time index to start
    alt_times = pre_spike(p).windows(1).all_windows;
    start_index = find(alt_times == start_sec);
    
    %% Get before manual rise times
    before_rise = pre_spike(p).windows(1).manual_before_rise;
    
    %% Define times to include
    % Going from start index to the last pre-rise time
    times_to_include = before_rise;
    if start_index > 1
        times_to_include(:,1:start_index-1) = 0;
    end
    
    %% Reduce size of times to include to align with the metric
    extra_size = length(alt_times) - length(times);
    times_to_include(:,1:extra_size) = [];
    test_alt_times = alt_times;
    test_alt_times(1:extra_size) = [];
    if ~isequal(round(test_alt_times*10)/10,round(times*10)/10)
        error('what'); 
    end
    if size(times_to_include,2) ~= length(times), error('what'); end
    
    %% Take the mode across spikes to get times for not spikes
    times_not_spikes = mode(times_to_include,1);
    times_not_spikes = repmat(times_not_spikes,...
        size(metrics.time.freq(1).(met).pt(p).not.data,1),1);
    
    %% Output into a struct
    include_times(p).times_to_include = logical(times_to_include);
    include_times(p).times_to_include_non_spike = logical(times_not_spikes);
    
end


end