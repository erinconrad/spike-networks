function network_change_fig

%{
I am not sure how to plot the NBS statistics, not sure what to show other
than a p-value which is kind of lame. Thinking I will just show the
permanova statistic.
%}

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
perm_folder = [results_folder,'perm_stats/'];


if exist(out_folder,'dir') == 0
    mkdir(out_folder);
end

% Load spike file for one patient to get the surround time
spike = load([eeg_folder,'HUP074_eeg.mat']);
spike = spike.spike;
surround_time = spike(1).surround_time;


% Loop through network types
listing = dir(perm_folder);
network_count = 0;
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
    
    network_folder = [perm_folder,name,'/'];
    
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
        sim = load([time_folder,pt_listing(1).name]);
        sim = sim.sim;
        nfreq = length(sim);
        
        stats(network_count).time(time_count).freq(1:nfreq).F_all = ...
            zeros(length(pt_listing),surround_time*2/time_window);
        
        stats(network_count).time(time_count).freq(1:nfreq).p_all = ...
            zeros(length(pt_listing),surround_time*2/time_window);
        
        % loop through pts
        for i = 1:length(pt_listing)
            
            pt_name = pt_listing(i).name;
            pt_name_pt = strsplit(pt_name,'_');
            pt_name_pt = pt_name_pt{1};
            
            % load pt file
            sim = load([time_folder,pt_name]);
            sim = sim.sim;
            
            for f = 1:nfreq
                stats(network_count).time(time_count).freq(f) = sim(f).name;
                stats(network_count).time(time_count).freq(f).F_all(i,:) = F;
                stats(network_count).time(time_count).freq(f).p_all(i,:) = p;
            end
            
            
        end
        
    end
        
    
    
end

end