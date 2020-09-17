function stats = get_perm

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

%% Now get F statistics for network differences
network_count = 0;
n_freq_abs = 0;
listing = dir(perm_folder);
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
        sim = load([time_folder,pt_listing(2).name]);
        sim = sim.sim;
        nfreq = length(sim);
        if n_freq_abs < nfreq
            n_freq_abs = nfreq;
        end
        
        

        all_names = {};
        % loop through pts
        for i = 1:length(pt_listing)

            fname = pt_listing(i).name;
            pt_name = strsplit(fname,'_');
            pt_name = pt_name{1};
            

            [a,b] = ismember(pt_name,all_names);
            if a == 1
                pt_idx = b;
            else
                all_names = [all_names;pt_name];
                pt_idx = length(all_names);
            end

            % load pt file
            sim = load([time_folder,fname]);
            sim = sim.sim;
            
            
            if isfield(sim(1),'index_windows') == 1 && isempty(sim(1).index_windows) == 0
                stats(network_count).time(time_count).index_windows = sim.index_windows;
                stats(network_count).time(time_count).fs = sim.fs;
            else
                fprintf('Warning, no index windows field for %s time window %1.1f. Skipping...\n',...
                    pt_name,time_window);
                continue;
            end
            %}

            
            
            % Add times
            temp_times = sim(1).index_windows(:,1)/sim(1).fs-3;
            times = realign_times(temp_times,nan);
            times = round(times*1e2)/(1e2); % rounding error
            stats(network_count).time(time_count).times = times;
            
            for f = 1:nfreq
                stats(network_count).time(time_count).freq(f).name = sim(f).name;

                
                % spike vs not a spike
                if contains(fname,'not') == 0
                    
                    if isfield(stats(network_count).time(time_count).freq(f),'F_all') && ...
                            size(stats(network_count).time(time_count).freq(f).F_all,2) > length(sim(f).F)
                        sim(f).F = [sim(f).F; ...
                            nan(size(stats(network_count).time(time_count).freq(f).F_all,2)-length(sim(f).F),1)];
                        fprintf('\nWarning, had to pad F stats with nans for %s.\n',pt_name);
                        % this is because I originally calculated adjacency
                        % matrices for -3 to +3 for some patients, but I do
                        % later patients and subsequent analyses on just -3
                        % to 0
                    end
                    stats(network_count).time(time_count).freq(f).F_all(pt_idx,:,1) = sim(f).F;
                else
                    stats(network_count).time(time_count).freq(f).F_all(pt_idx,:,2) = sim(f).F;
                end
                
                % These pseudo F stats are from a permANOVA comparing the
                % network in the first time period to that of subsequent
                % time periods
                
                
            end


        end

    end

end


end