function stats = get_all_metrics(windows)

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
sp_diff_folder = [results_folder,'net_diff_stats/'];

%% Now get F statistics for network differences
network_count = 0;
n_freq_abs = 0;
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
        
        % Skip if I didn't ask for it
        if ~ismember(time_window,windows), continue; end

        time_count = time_count + 1;
        stats(network_count).time(time_count).name = time_name;
        stats(network_count).time(time_count).time_window = time_window;
        
        time_folder = [network_folder,time_name,'/'];

        pt_listing = dir([time_folder,'*.mat']);

        % load one to get nfreq
        sim = load([time_folder,pt_listing(2).name]);
        sim = sim.metrics;
        nfreq = length(sim.freq);
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
            sim = sim.metrics;

            
            % Also load the spike diff file if it exists
            if contains(fname,'not')
                nstext = '_not_spike';
            else
                nstext = '';
            end
            sp_diff_file = [sp_diff_folder,...
                stats(network_count).name,'/',time_name,'/',...
                pt_name,nstext,'_perm.mat'];
            if exist(sp_diff_file,'file') ~= 0
                add_sp_diff = 1;
                spd = load(sp_diff_file);
                spd = spd;
            else
                add_sp_diff = 0;
            end
                
            
            for f = 1:nfreq
                stats(network_count).time(time_count).freq(f).name = sim.freq(f).name;
                ns = sim.freq(f).ns.data;
                ns_name = sim.freq(f).ns.name;
                ge = sim.freq(f).ge.data;
                ge_name = sim.freq(f).ge.name;
                ns_all = sim.freq(f).ns_all.data;

                
                % Avg ns_all across channels
                ns_avg = squeeze(mean(ns_all,3));
                ns_avg_name = 'average node strength';
                
                stats(network_count).time(time_count).freq(f).ge.pt(pt_idx).name = pt_name;
                stats(network_count).time(time_count).freq(f).ns_avg.pt(pt_idx).name = pt_name;
                stats(network_count).time(time_count).freq(f).ns_big.pt(pt_idx).name = pt_name;
                
                % spike vs not a spike
                if contains(fname,'not') == 1
                    stats(network_count).time(time_count).freq(f).ge.pt(pt_idx).not.data(:,:) = ge;
                    stats(network_count).time(time_count).freq(f).ns_avg.pt(pt_idx).not.data(:,:) = ns_avg;
                    stats(network_count).time(time_count).freq(f).ns_big.pt(pt_idx).not.data(:,:) = ns;
                else
                    stats(network_count).time(time_count).freq(f).ge.pt(pt_idx).spike.data(:,:) = ge;
                    stats(network_count).time(time_count).freq(f).ns_avg.pt(pt_idx).spike.data(:,:) = ns_avg;
                    stats(network_count).time(time_count).freq(f).ns_big.pt(pt_idx).spike.data(:,:) = ns;
                end
                
                stats(network_count).time(time_count).freq(f).ns_big.name = 'Node strength (spike channel)';
                stats(network_count).time(time_count).freq(f).ns_avg.name = 'Node strength (average)';
                stats(network_count).time(time_count).freq(f).ge.name = 'Global efficiency';
                
                % Add spike diff stuff
                if add_sp_diff == 1
                    stats(network_count).time(time_count).freq(f).F.index_windows = spd.sim(f).index_windows;
                    stats(network_count).time(time_count).freq(f).score.index_windows = spd.sim(f).index_windows;
                    % convert score to a 2 dimensional vector
                    score = spd.sim(f).score;
                    time_idx = spd.sim(f).time_idx;
                    num_time_idx = length(unique(time_idx));
                    new_score = zeros(length(score)/num_time_idx,num_time_idx);
                    for tt = 1:num_time_idx
                        new_score(:,tt) = score(time_idx == tt);
                    end
                    
                    if contains(fname,'not') == 1
                        stats(network_count).time(time_count).freq(f).F.pt(pt_idx).not.data = spd.sim(f).F';
                        stats(network_count).time(time_count).freq(f).score.pt(pt_idx).not.data = new_score;
                    else
                        stats(network_count).time(time_count).freq(f).F.pt(pt_idx).spike.data = spd.sim(f).F';
                        stats(network_count).time(time_count).freq(f).score.pt(pt_idx).spike.data = new_score;
                    end
                end
                
            end


        end

    end

end

end