function net_fig(windows)

%% Parameters
alpha = 0.05;

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


%% First, get signal deviations
% get full directory listing
listing = dir(sig_dev_folder);
count = 0;
for i = 1:length(listing)
    % look for only those that are directories, ignoring '.' and '..'
    if listing(i).isdir == 0
        continue;
    end
    
    if strcmp(listing(i).name,'.') == 1 || strcmp(listing(i).name,'..') == 1
        continue
    end
    
    
    
    % assume all other things are directories with different time windows
    
    time_text = listing(i).name;
    time_window = str2num(time_text);
    time_window_folder = [sig_dev_folder,time_text,'/'];
    
    % load the file
    sub_listing = dir(time_window_folder);
    for k = 1:length(sub_listing)
        
        if contains(sub_listing(k).name,'.mat') == 0, continue; end
        
        
        if contains(sub_listing(k).name,'not_spike') == 1, continue; end

        count = count+1;
        temp_sig_dev = load([time_window_folder,sub_listing(k).name]);
        sig_dev(count).name = time_text;
        sig_dev(count).time_window = time_window;
        sig_dev(count).sig_dev = temp_sig_dev.sig_dev;
    end
end

%% Now combine t statistics across patients to find the times of significant power change
n_windows = count;
for t = 1:n_windows
    sig_dev(t).t_stat_all = nan(length(sig_dev(1).sig_dev),...
        length(sig_dev(1).sig_dev(1).stats));
    for i = 1:length(sig_dev(t).sig_dev)
        for tt = 1:length(sig_dev(t).sig_dev(i).stats)
            if isempty(sig_dev(t).sig_dev(i).stats(tt).tstat) == 0
                sig_dev(t).t_stat_all(i,tt) = sig_dev(t).sig_dev(i).stats(tt).tstat;
            end
        end
    end
    
    sig_dev(t).p_all = nan(size(sig_dev(t).t_stat_all,2),1);
    sig_dev(t).sig = zeros(size(sig_dev(t).t_stat_all,2),1);
    for tt = 2:size(sig_dev(t).t_stat_all,2)
        [~,p] = ttest(sig_dev(t).t_stat_all(:,tt));
        sig_dev(t).p_all(tt) = p;
        if p < alpha/(size(sig_dev(t).t_stat_all,2)-1)
            sig_dev(t).sig(tt:end) = 1;
            break
        end
    end
    sig_dev(t).sig(1) = 1; % to ignore the first time
    sig_dev(t).sig = logical(sig_dev(t).sig);
end

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
            
            if isfield(sim(1),'index_windows') == 1 && isempty(sim(1).index_windows) == 0
                stats(network_count).time(time_count).index_windows = sim.index_windows;
                stats(network_count).time(time_count).fs = sim.fs;
            end


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

%% Now plot network change over time and look for significant slope
figure
set(gcf,'position',[1 100 1399 600])
[ha, ~] = tight_subplot(time_count, n_freq_abs+1, [0.1 0.01], [0.12 0.07], [0.05 0.01]);
z_range = zeros(time_count,2);
for n = 1:network_count

    net_name = stats(n).name;
    
    for t = 1:time_count
        
        if t>size(stats(n).time), continue; end
        
        % change times for x axis
        nchunks = size(stats(n).time(t).freq(1).F_all,2);
        
        if isfield(stats(n).time(t),'index_windows') && isempty(stats(n).time(t).index_windows) == 0
            temp_times = stats(n).time(t).index_windows(:,1)/stats(n).time(t).fs-3;
            times = realign_times(temp_times,surround_time);
        else

            times = realign_times(nchunks,surround_time);
        end
        sig_power_change_times = sig_dev(t).sig;
        times = times(~sig_power_change_times);
        
        nfreq = length(stats(n).time(t).freq);
        for f = 1:nfreq
            
            % Get appropriate subplot
            if strcmp(net_name,'coherence') == 1
                column_add = 1;
            else
                column_add = 0;
            end
            % this adds the number of frequencies + 1 if it's on the 2nd
            % time point (to move down a row), and it adds which frequency
            % (which is 1 if simple) and adds 1 if coherence, to start with
            % the 2nd column for coherence
            %sp = f + column_add;
            sp = (n_freq_abs+1)*(t-1) + f + column_add;
            axes(ha(sp));
            
            % Get F stats
            F = stats(n).time(t).freq(f).F_all;
            
            % Just take times without significant power change
            F_curr = F(:,~sig_power_change_times);
            
            
            % z score to normalize within pt
            z_curr = (F_curr-nanmean(F_curr,2))./nanstd(F_curr,0,2); %nan because first time is nans
            slopes = zeros(size(z_curr,1),1);
            
            % loop over patients and plot
            for i = 1:size(z_curr,1)
                plot(times,z_curr(i,:),'ko');
                hold on
                
                % Get slope for each patient
                y = z_curr(i,:)';
                x = [ones(length(y),1), (1:length(y))'];
                % do regression to find best fit line through the F stats
                % for that patient
                b = x\y; 
                slopes(i) = b(2);
                
                
            end
            
            % See if slopes are significantly different from zero
            [~,p] = ttest(slopes);
            
            % Plot an overall trend line
           % x = repmat(1:size(z_curr,2),size(z_curr,1),1);
            x = repmat(times',size(z_curr,1),1);
            y = z_curr(:);
            x = x(:);
            X = [ones(length(x),1),x];
            b = X\y;
            
            if p < alpha/(n_freq_abs+1)
                %plot(times',b(1)+b(2)*(1:size(z_curr,2)),'g','linewidth',2);
                plot(times',b(1)+b(2)*times,'g','linewidth',2);
            else
                %plot(times',b(1)+b(2)*(1:size(z_curr,2)),'k','linewidth',2);
                plot(times',b(1)+b(2)*times,'k','linewidth',2);
            end
            
            set(gca,'fontsize',20)
            if f == 4 && t == time_count
                xlabel('Time relative to spike peak (s)')
            end 
            if n == 2
                ylabel(sprintf('Network distance from\nfirst time (z-score)'))
            end
            
            % adjust z_range
            if max(max(z_curr)) > z_range(t,2)
                z_range(t,2) = max(max(z_curr));
            end

            if min(min(z_curr)) < z_range(t,1)
                z_range(t,1) = min(min(z_curr));
            end
            
        end
        
    end
    
end
for sp = 1:length(ha)
    axes(ha(sp))
    t = ceil(sp/(n_freq_abs+1));
    ylim(z_range(t,:))
    
    if mod(sp,(n_freq_abs+1)) ~= 1
        yticklabels([])
    end
end

%% Now, do the same thing for the power change (should not have a slope)
figure
set(gcf,'position',[1 100 1399 600])
[ha, ~] = tight_subplot(time_count, 1, [0.1 0.01], [0.12 0.07], [0.05 0.01]);
n = 1;
for t = 1:time_count
    % change times for x axis
    nchunks = size(stats(n).time(t).freq(1).F_all,2);

    if isfield(stats(n).time(t),'index_windows') && isempty(stats(n).time(t).index_windows) == 0
        temp_times = stats(n).time(t).index_windows(:,1)/stats(n).time(t).fs-3;
        times = realign_times(temp_times,surround_time);
    else

        times = realign_times(nchunks,surround_time);
    end
    sig_power_change_times = sig_dev(t).sig;
    times = times(~sig_power_change_times);
    
    axes(ha(t));
    
    % get t-stats
    t_stats = sig_dev(t).t_stat_all;
    
    % Just take times without significant power change           
    t_curr = t_stats(:,~sig_power_change_times);
    
    % z score to normalize within pt
    z_curr = (t_curr-nanmean(t_curr,2))./nanstd(t_curr,0,2); %nan because first time is nans
    slopes = zeros(size(t_curr,1),1);
    
    % loop over patients and plot
    for i = 1:size(z_curr,1)
        plot(times,z_curr(i,:),'ko');
        hold on

        % Get slope for each patient
        y = z_curr(i,:)';
        x = [ones(length(y),1), (1:length(y))'];
        % do regression to find best fit line through the F stats
        % for that patient
        b = x\y; 
        slopes(i) = b(2);
    end

    % See if slopes are significantly different from zero
    [~,p] = ttest(slopes);

    % Plot an overall trend line
    %x = repmat(1:size(z_curr,2),size(z_curr,1),1);
    x = repmat(times',size(z_curr,1),1);
    y = z_curr(:);
    x = x(:);
    x = [ones(length(x),1),x];
    b = x\y;

    if p < alpha
        plot(times,b(1)+b(2)*times,'g','linewidth',2);
    else
        plot(times,b(1)+b(2)*times,'k','linewidth',2);
    end
    
    set(gca,'fontsize',20)
    
    xlabel('Time relative to spike peak (s)')
    
    ylabel(sprintf('Network distance from\nfirst time (z-score)'))
    
end

end