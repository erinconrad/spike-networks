function pre_spike = find_pre_spike_windows(times)

%% Get file locations, load spike times and pt structure
locations = spike_network_files;
main_folder = locations.main_folder;
results_folder = [main_folder,'results/'];
outer_spike_rise_folder = [results_folder,'two_reviewers_rise/'];

%% Load spike rise structs
for r = 1:2
    spike_rise_folder = [outer_spike_rise_folder,'Rev',num2str(r),'/'];
    listing = dir(spike_rise_folder);
    pt_idx = 0;
    for i = 1:length(listing)

        fname = listing(i).name;

        % Skip if . or ..
        if strcmp(fname,'.') == 1 || strcmp(fname,'..') == 1 || strcmp(fname,'.DS_Store') == 1
            continue
        end
        pt_idx = pt_idx + 1;

        early = load([spike_rise_folder,fname]);
        early = early.early;
        name = early.name;
        pre_spike(pt_idx).name = name;

        % Initialize if it's the first reviewer
        if r == 1
            pre_spike(pt_idx).times = times';
            
            pre_spike(pt_idx).before_rise = ...
                zeros(length(early.spike),length(pre_spike(pt_idx).times));
            pre_spike(pt_idx).before_rise_both = ...
                zeros(length(early.spike),length(pre_spike(pt_idx).times),2);
            pre_spike(pt_idx).before_rise_both_times = ...
                zeros(length(early.spike),2);
        end

        % Loop through spikes
        for s = 1:length(early.spike)

            % Get the time (relative to peak) of the earliest change
            change_time = early.spike(s).time;
            
            % Fix for bug in getting Jim's times
            %{
            I didn't reinitialize the structure array for each patient and
            so later patients with fewer spikes will have extra spikes at
            the end filled in from prior patients with more spikes. I will
            ignore any spikes beyond the size of Erin's array.
            %}
            if r == 2 && s > size(pre_spike(pt_idx).before_rise,1)
                break
            end
                
            
            pre_spike(pt_idx).before_rise_both_times(s,r) = change_time;

            
            time_window = times(2) - times(1);

            % I want the end of the window to be earlier than the time of 
            % earliest spike change
            window_end = times' + time_window; 

            % Find windows in which the end time is before the time of
            % earliest spike change
            windows_before_rise = (window_end < change_time)';

            % Store this logical index
            pre_spike(pt_idx).before_rise_both(s,:,r) = windows_before_rise;

            % If r == 2, pull up reviewer 1
            if r == 2
                windows_before_rise_1 = pre_spike(pt_idx).before_rise_both(s,:,1);



                % Find windows that occur between rise for both reviewers
                windows_before_rise_new = windows_before_rise & windows_before_rise_1;
                % So, if either one says the window is not before the
                % rise (windows_before_rise == 0) then this will also
                % say it's not before the rise (be 0)

                pre_spike(pt_idx).before_rise(s,:) = windows_before_rise_new;
            end
                

                

           

        end
    end
end

end