function eval_slope(sd,not_a_spike,windows)

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
nbs_folder = [results_folder,'nbs_stats/'];

freq_names = {'delta','theta','alpha','beta','low_gamma',...
    'high_gamma','ultra_high','broadband'};

F_times = [0.1 2 8;...
    0.2 2 9;...
    0.5 2 5;...
    1 2 2];

if exist(out_folder,'dir') == 0
    mkdir(out_folder);
end

% Load spike file for one patient to get the surround time
spike = load([eeg_folder,'HUP074_eeg.mat']);
spike = spike.spike;
surround_time = spike(1).surround_time;


% Loop through network types
if sd == 0
    network_count = 0;
    n_freq_abs = 0;
    max_F = 0;
    max_z = 0;
    min_z = 0;
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
            
            % Skip if not the time window we want
            if ismember(time_window,windows) == 0, continue; end

            time_count = time_count + 1;
            stats(network_count).time(time_count).name = time_name;
            stats(network_count).time(time_count).time_window = time_window;
            time_folder = [network_folder,time_name,'/'];

            pt_listing = dir([time_folder,'*.mat']);

            % load one to get nfreq
            sim = load([time_folder,pt_listing(1).name]);
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


                for f = 1:nfreq
                    stats(network_count).time(time_count).freq(f).name = sim(f).name;
                    
                    if contains(fname,'not') == 0
                        stats(network_count).time(time_count).freq(f).F_all(pt_idx,:,1) = sim(f).F;
                    else
                    	stats(network_count).time(time_count).freq(f).F_all(pt_idx,:,2) = sim(f).F;
                    end
                end


            end

        end

    end
else
    n_freq_abs = 0;

    % Loop through time scales
    time_listing = dir(sig_dev_folder);
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
        
        % Skip if not the time window we want
        if ismember(time_window,windows) == 0, continue; end

        time_count = time_count + 1;
        stats.time(time_count).name = time_name;
        stats.time(time_count).time_window = time_window;
        time_folder = [sig_dev_folder,time_name,'/'];

        pt_listing = dir([time_folder,'*.mat']);
        
        for i = 1:length(pt_listing)            
            % load pt file
            sig_dev = load([time_folder,pt_listing(i).name]);
            sig_dev = sig_dev.sig_dev;
            
            for p = 1:length(sig_dev)
                if contains(pt_listing(i).name,'not') == 1
                    stats.time(time_count).freq.z_all(p,:,2) = sig_dev(p).z_score_dev;
                else
                    stats.time(time_count).freq.z_all(p,:,1) = sig_dev(p).z_score_dev;
                end
            end
            
        end
        
        
    end
    
end


for n = 1:length(stats)
    for t = 1:length(stats(n).time)
        time_window = stats(n).time(time_count).time_window;
        for f = 1:length(stats(n).time(t).freq)
            if sd == 1
                F_all = stats(n).time(t).freq(f).z_all;
            else
                F_all = stats(n).time(t).freq(f).F_all;
            end
            
            if not_a_spike == 1
                F_all = F_all(:,:,2);
            else
                F_all = F_all(:,:,1);
            end
            
            n_patients = size(F_all,1);
            slopes = zeros(n_patients,1);
            
            % restrict times to only those pre-spike
            [~,b] = ismember(time_window,F_times);
            F_res = squeeze(F_all(:,F_times(b,2):F_times(b,end),1));
            
            for i = 1:n_patients
                
                
                y = F_res(i,:)';
                y = (y-nanmean(y))./nanstd(y);
                x = [ones(length(y),1), (1:length(y))'];
                b = x\y;
                
                % do regression to find best fit line through the F stats
                % for that patient
                slopes(i) = b(2);
                
            end
            
            all_slopes(n).time(t).freq(f).slope = slopes;
            [~,p] = ttest(slopes);
            all_slopes(n).time(t).freq(f).p = p;
        end
    end
end

end