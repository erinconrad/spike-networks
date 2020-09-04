function slope_pca(sd,windows)

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


%% Do PCA on F statistic data
for t = 1:time_count
    figure
    set(gcf,'position',[1 100 1399 500])
    [ha, pos] = tight_subplot(2, n_freq_abs+1, [0.01 0.01], [0.12 0.07], [0.05 0.01]);
    for n = 1:length(stats)

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
            
            for s = 1:2
                % this adds the number of frequencies + 1 if it's on the 2nd
                % time point (to move down a row), and it adds which frequency
                % (which is 1 if simple) and adds 1 if coherence, to start with
                % the 2nd column for coherence
                sp = (n_freq_abs+1)*(s-1) + f + column_add;
                axes(ha(sp));

                
                if sd == 0
                    if ndims(stats(n).time(t).freq(f).F_all) < s+1, continue; end
                    F = squeeze(stats(n).time(t).freq(f).F_all(:,:,s)); % n_patients x n_times
                else
                    if ndims(stats(n).time(t).freq(f).z_all) < s+1, continue; end
                    F = squeeze(stats(n).time(t).freq(f).z_all(:,:,s));
                end

                [a,b] = ismember(stats(n).time(time_count).time_window,F_times);
                
                
                % only take times 1:before signficant signal deviation
                F = F(:,F_times(b,2):F_times(b,end));

                % Take a look at correlation coefficient
                R = corrcoef(F);

                % Do PCA
                [coeff,score,latent] = pca(F);

                % Plot the coefficients of the first component
                plot(coeff(:,1))
                hold on

                
                if s == 2 && f == 4
                     xlabel('Time period')
                end 

                if s == 2 && f == 1 && n == 2
                    ylabel('Coefficient of first component')
                end
                %}

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

%% Test if the scores for spike are higher than the scores for not a spike
for t = 1:time_count
    for n = 1:length(stats)
        for f = 1:length(stats(n).time(t).freq)
           
            if sd == 0
                F = (stats(n).time(t).freq(f).F_all); %n_patients x n_time periods x 2 (spike and not a spike)
            else
                F = (stats(n).time(t).freq(f).z_all);
            end
            
            % Concatenate spike and not a spike
            is_spike_idx = [ones(size(F,1),1);zeros(size(F,1),1)];       
            F_cat = [squeeze(F(:,:,1));squeeze(F(:,:,2))];
            
            
            
        end
    end
end



end