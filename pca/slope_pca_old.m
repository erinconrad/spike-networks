function slope_pca(sd,not_a_spike)

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
            
            % Skip if not 0.5 or 0.2
            if contains(time_name,'0.5') == 0 && contains(time_name,'0.2')== 0
                continue
            end

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


           if isfield(sim,'index_windows')



               for f = 1:nfreq
                    stats(network_count).time(time_count).index_windows = sim(f).index_windows;
                    stats(network_count).time(time_count).fs = sim(f).fs;

                    stats(network_count).time(time_count).freq(f).F_all = ...
                        nan(length(pt_listing),size(sim(f).index_windows,1));

                    stats(network_count).time(time_count).freq(f).p_all = ...
                        nan(length(pt_listing),size(sim(f).index_windows,1));

                    stats(network_count).time(time_count).freq(f).z_all = ...
                        nan(length(pt_listing),size(sim(f).index_windows,1));

                end
           else

               for f = 1:nfreq
                    stats(network_count).time(time_count).freq(f).F_all = ...
                        nan(length(pt_listing),surround_time*2/time_window);

                    stats(network_count).time(time_count).freq(f).p_all = ...
                        nan(length(pt_listing),surround_time*2/time_window);

                    stats(network_count).time(time_count).freq(f).z_all = ...
                        nan(length(pt_listing),surround_time*2/time_window);

                end
           end



            pt_names = {};

            % loop through pts
            for i = 1:length(pt_listing)

                pt_name = pt_listing(i).name;
                pt_name_pt = strsplit(pt_name,'_');
                pt_name_pt = pt_name_pt{1};
                pt_names = [pt_names;pt_name_pt];

                % load pt file
                sim = load([time_folder,pt_name]);
                sim = sim.sim;


                for f = 1:nfreq
                    stats(network_count).time(time_count).freq(f).name = sim(f).name;
                    stats(network_count).time(time_count).freq(f).F_all(i,:) = sim(f).F;
                    stats(network_count).time(time_count).freq(f).p_all(i,:) = sim(f).p;

                    F_curr = stats(network_count).time(time_count).freq(f).F_all;
                    if max_F < max(max(F_curr))
                        max_F = max(max(F_curr));
                    end



                    % convert F stats to z scores to compare across time points
                    stats(network_count).time(time_count).freq(f).z_all(i,:) = (sim(f).F-nanmean(sim(f).F))./nanstd(sim(f).F);

                    if max_z < max(max(stats(network_count).time(time_count).freq(f).z_all))
                        max_z = max(max(stats(network_count).time(time_count).freq(f).z_all));
                    end

                    if min_z > min(min(stats(network_count).time(time_count).freq(f).z_all))
                        min_z = min(min(stats(network_count).time(time_count).freq(f).z_all));
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
        
        % Skip if not 0.5 or 0.2
            if contains(time_name,'0.5') == 0 && contains(time_name,'0.2')== 0
                continue
            end

        time_count = time_count + 1;
        stats.time(time_count).name = time_name;
        stats.time(time_count).time_window = time_window;
        time_folder = [sig_dev_folder,time_name,'/'];

        pt_listing = dir([time_folder,'*.mat']);
        
        for i = 1:length(pt_listing)
            if contains(pt_listing(i).name,'not_spike') == 1
                if not_a_spike == 0, continue; end
            else
                if not_a_spike == 1, continue; end
            end
            
            % load pt file
            sig_dev = load([time_folder,pt_listing(i).name]);
            sig_dev = sig_dev.sig_dev;
            
            for p = 1:length(sig_dev)
                stats.time(time_count).freq.z_all(p,:) = sig_dev(i).z_score_dev;
            end
            
        end
        
        
    end
    
end


%% Do PCA on F statistic data
figure
set(gcf,'position',[1 100 1399 500])
[ha, pos] = tight_subplot(time_count, n_freq_abs+1, [0.01 0.01], [0.12 0.07], [0.05 0.01]);
for n = 1:length(stats)
    for t = 1:time_count
        if t == 1
            F_times = [2:8];
        elseif t == 2
            F_times = [2:5];
        end
        
        for f = 1:length(stats(n).time(t).freq)
            
            % Get appropriate subplot
            if sd == 0
            if strcmp(stats(n).name,'coherence') == 1
                column_add = 1;
            else
                column_add = 0;
            end
            else
                column_add = 0;
            end
            % this adds the number of frequencies + 1 if it's on the 2nd
            % time point (to move down a row), and it adds which frequency
            % (which is 1 if simple) and adds 1 if coherence, to start with
            % the 2nd column for coherence
            sp = (n_freq_abs+1)*(t-1) + f + column_add;
            axes(ha(sp));
            
            if sd == 0
                F = stats(n).time(t).freq(f).F_all; % n_patients x n_times
            else
                F = stats(n).time(t).freq(f).z_all;
            end
            
            % only take times 2:before signficant signal deviation
            F = F(:,F_times);
            
            % Take a look at correlation coefficient
            %R = corrcoef(F);
            
            % Do PCA
            [coeff,score,latent] = pca(F);
            
            % Plot the coefficients of the first component
            plot(coeff(:,1))
            hold on
            
            if t == 3 && f == 4
                 xlabel('Time period')
            end 
            
            if t == 2 && f == 1 && n == 2
                ylabel('Coefficient of first component')
            end
            
            xticklabels([])
            yticklabels([])
            
            if sd == 0
            if t == 1 && strcmp(stats(n).name,'coherence') == 1
                title(sprintf('%s',...
                    strrep(stats(n).time(t).freq(f).name,'_',' ')))
            elseif t == 1 && strcmp(stats(n).name,'simple') == 1
                title('correlation')
            end
            end

        end
    end
end


end