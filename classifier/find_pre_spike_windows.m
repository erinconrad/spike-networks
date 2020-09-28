function pre_spike = find_pre_spike_windows(windows)

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

%% Load spike rise structs
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

    
    for t = 1:length(windows)
        pre_spike(pt_idx).windows(t).which = windows(t);
        if windows(t) == 0.1
            pre_spike(pt_idx).windows(t).all_windows = [(-1:windows(t):0)'];
        elseif windows(t) == 0.2
            pre_spike(pt_idx).windows(t).all_windows = [(-2:windows(t):0)'];
        elseif windows(t) == 0.5
            pre_spike(pt_idx).windows(t).all_windows = [(-3:windows(t):0)'];
        end
        pre_spike(pt_idx).windows(t).before_rise = ...
            zeros(length(early.spike),length(pre_spike(pt_idx).windows(t).all_windows));
    end
    
    % Loop through spikes
    for s = 1:length(early.spike)
        
        % Get the time (relative to peak) of the earliest change
        change_time = early.spike(s).time;
        
        % Loop through windows
        for t = 1:length(windows)
            all_windows = pre_spike(pt_idx).windows(t).all_windows;
            
            % I want the end of the window to be earlier than the time of 
            % earliest spike change
            window_end = all_windows + windows(t); 
            
            % Find windows in which the end time is before the time of
            % earliest spike change
            windows_before_rise = window_end < change_time;
            
            % Store this logical index
            pre_spike(pt_idx).windows(t).before_rise(s,:) = windows_before_rise';
            
        end
        
    end
end

end