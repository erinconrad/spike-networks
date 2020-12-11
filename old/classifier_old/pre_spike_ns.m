function nets = pre_spike_ns(pre_spike,windows,met)

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

% Loop through network types
network_count = 0;
listing = dir(ns_folder);
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

    network_folder = [ns_folder,name,'/'];

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
        
        % Skip if I didn't request it
        if ~ismember(time_window,windows), continue; end

        time_count = time_count + 1;
        stats(network_count).time(time_count).name = time_name;
        stats(network_count).time(time_count).time_window = time_window;
        
        time_folder = [network_folder,time_name,'/'];

        pt_listing = dir([time_folder,'*.mat']);
        
        all_names = {};
        % Loop through pts
        for i = 1:length(pt_listing)

            fname = pt_listing(i).name;
            pt_name = strsplit(fname,'_');
            pt_name = pt_name{1};
            
            % Skip it if it doesn't have a corresponding pre-spike entry
            sp_idx = 0;
            for j = 1:length(pre_spike)
                sp_name = pre_spike(j).name;
                if strcmp(sp_name,pt_name)
                    sp_idx = j;
                    break
                end
            end
            if sp_idx == 0, continue; end
            
            
            [a,b] = ismember(pt_name,all_names);
            if a == 1
                pt_idx = b;
            else
                all_names = [all_names;pt_name];
                pt_idx = length(all_names);
            end
        
            % Load pt file
            metrics = load([time_folder,fname]);
            metrics = metrics.metrics;
            
            % Get the corresponding time windows to keep from the
            % pre_spike structure
            before_rise = pre_spike(sp_idx).windows(time_count).before_rise;
            
            % Loop through spikes
            %for s = 1:length(meta.spike)
                
                
                
                
            %end
            
        end
        
    end
    
end



end