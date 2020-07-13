function plot_metrics(simple,time_window)

%% Get file locations, load spike times and pt structure
locations = spike_network_files;
main_folder = locations.main_folder;
results_folder = [main_folder,'results/'];
data_folder = [main_folder,'data/'];
eeg_folder = [main_folder,'results/eeg_data/'];
script_folder = locations.script_folder;
addpath(genpath(script_folder));
pt_file = [data_folder,'spike_structures/pt.mat'];
bct_folder = locations.BCT;
addpath(genpath(bct_folder));
out_folder = [results_folder,'plots/'];
sig_dev_folder = [results_folder,'signal_deviation/manual/'];
metrics_folder = [results_folder,'metrics/'];


if exist(out_folder,'dir') == 0
    mkdir(out_folder);
end

% Load spike file for one patient to get the surround time
spike = load([eeg_folder,'HUP074_eeg.mat']);
spike = spike.spike;
surround_time = spike(1).surround_time;

% Loop through network types
listing = dir(metrics_folder);
network_count = 0;
n_freq_abs = 0;
max_F = 0;
for l = 1:length(listing)
    name= listing(l).name;
    
    % Skip if . or ..
    if strcmp(name,'.') == 1 || strcmp(name,'..') == 1
        continue
    end
    
    % Skip if not a directory
    if listing(l).isdir == 0, continue; end
    
    network_count = network_count + 1;
    stats(network_count).name = name;
    
    network_folder = [metrics_folder,name,'/'];
    
    % Loop through time scales
    time_listing = dir(network_folder);
    time_count = 0;
    
    for k = 1:length(time_listing)
        time_name= time_listing(k).name;
        time_window = str2num(time_name);
        
        % Skip if . or ..
        if strcmp(time_name,'.') == 1 || strcmp(time_name,'..') == 1
            continue
        end

        % Skip if not a directory
        if time_listing(k).isdir == 0, continue; end
        
        time_count = time_count + 1;
        stats(network_count).time(time_count).name = time_name;
        stats(network_count).time(time_count).time_window = time_window;
        time_folder = [network_folder,time_name,'/'];
        
        pt_listing = dir([time_folder,'*.mat']);
        
        % load one to get nfreq
        metrics = load([time_folder,pt_listing(1).name]);
        metrics = metrics.metrics;
        nfreq = length(metrics);
        if n_freq_abs < nfreq
            n_freq_abs = nfreq;
        end
        
        for f = 1:nfreq
            stats(network_count).time(time_count).freq(f).ge.data = ...
                zeros(length(pt_listing),surround_time*2/time_window);

            stats(network_count).time(time_count).freq(f).p_all = ...
                zeros(length(pt_listing),surround_time*2/time_window,2);
        end
        
        % loop through pts
        for i = 1:length(pt_listing)
            
            pt_name = pt_listing(i).name;
            pt_name_pt = strsplit(pt_name,'_');
            pt_name_pt = pt_name_pt{1};
            
            % load pt file
            metrics = load([time_folder,pt_name]);
            metrics = metrics.metrics;
            involved = metrics.involved;
            
            for f = 1:nfreq
                stats(network_count).time(time_count).freq(f).name = metrics.freq(f).name;
                stats(network_count).time(time_count).freq(f).ge.name = metrics.freq(f).ge.name;
                stats(network_count).time(time_count).freq(f).ge.data(i,:) = metrics.freq(f).ge.data;
                stats(network_count).time(time_count).freq(f).ns.name = metrics.freq(f).ns.name;
                stats(network_count).time(time_count).freq(f).bc.name = metrics.freq(f).bc.name;
                
                % Get ns, bc for involved and uninvolved chs
                stats(network_count).time(time_count).freq(f).ns.data(i,:,1) = metrics.freq(f).ns.data;
                
                stats(network_count).time(time_count).freq(f).bc.data(i,:,2) = metrics.freq(f).bc.data;
                
                F_curr = stats(network_count).time(time_count).freq(f).F_all;
                if max_F < max(max(F_curr))
                    max_F = max(max(F_curr));
                end
            
            end
            
            
        end
        
    end
 
end


end