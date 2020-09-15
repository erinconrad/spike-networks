function ns_fig(windows)

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
ns_folder = [results_folder,'metrics/manual/'];

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
listing = dir(ns_folder);
network_count = 0;
n_freq_abs = 0;
max_z = 0;
min_z = 0;
network_names = {};
for l = 1:length(listing)
    name= listing(l).name;
    
    % Skip if . or ..
    if strcmp(name,'.') == 1 || strcmp(name,'..') == 1|| strcmp(name,'.DS_Store') == 1
        continue
    end
    
    % Skip if not a directory
    if listing(l).isdir == 0, continue; end
    
    network_names = [network_names;name];
    
    network_count = network_count + 1;
    stats(network_count).name = name;
    
    network_folder = [ns_folder,name,'/'];
    
    % Loop through time scales
    time_listing = dir(network_folder);
    time_count = 0;
    time_names = {};
    
    for k = 1:length(time_listing)
        time_name= time_listing(k).name;
        time_window = str2num(time_name);
        
        % Skip if . or ..
        if strcmp(time_name,'.') == 1 || strcmp(time_name,'..') == 1|| strcmp(name,'.DS_Store') == 1
            continue
        end

        % Skip if not a directory
        if time_listing(k).isdir == 0, continue; end
        
        % Skip if not the time window we want
        if ismember(time_window,windows) == 0, continue; end

        
        time_names = [time_names;time_name];
        
        time_count = time_count + 1;
        stats(network_count).time(time_count).name = time_name;
        stats(network_count).time(time_count).time_window = time_window;
        time_folder = [network_folder,time_name,'/'];
        
        pt_listing = dir([time_folder,'*.mat']);
        
        % load one to get nfreq and index windows
        sim = load([time_folder,pt_listing(1).name]);
        sim = sim.metrics;
        nfreq = length(sim.freq);
        if n_freq_abs < nfreq
            n_freq_abs = nfreq;
        end
        stats(network_count).time(time_count).index_windows = sim.index_windows;
        stats(network_count).time(time_count).fs = sim.fs;
       
        
        for f = 1:nfreq
            stats(network_count).time(time_count).freq(f).ns_all = ...
                nan(length(pt_listing),size(sim.index_windows,1));
            stats(network_count).time(time_count).freq(f).z_all = ...
                nan(length(pt_listing),size(sim.index_windows,1));
            stats(network_count).time(time_count).freq(f).p_all = ...
                nan(length(pt_listing),size(sim.index_windows,1));
            stats(network_count).time(time_count).freq(f).t_all = ...
                nan(length(pt_listing),size(sim.index_windows,1));
        end
        
        
        pt_names = {};
        % loop through pts
        for i = 1:length(pt_listing)
            
            pt_name = pt_listing(i).name;
            
            if contains(pt_name,'ns') == 0, continue; end
            
            pt_name_pt = strsplit(pt_name,'_');
            pt_name_pt = pt_name_pt{1};
            pt_names = [pt_names;pt_name_pt];
            
            % load pt file
            sim = load([time_folder,pt_name]);
            sim = sim.metrics;
            
            
            for f = 1:nfreq
                stats(network_count).time(time_count).freq(f).name = sim.freq(f).name;
                ns_all = sim.freq(f).ns.data;
                ns_avg = squeeze(mean(ns_all,1));
                stats(network_count).time(time_count).freq(f).ns_all(i,:) = ns_avg;
              
                % convert node strengths to z scores to compare across time points
                stats(network_count).time(time_count).freq(f).z_all(i,:) = (ns_avg-nanmean(ns_avg))./nanstd(ns_avg);
            
                if max_z < max(max(stats(network_count).time(time_count).freq(f).z_all))
                    max_z = max(max(stats(network_count).time(time_count).freq(f).z_all));
                end
                
                if min_z > min(min(stats(network_count).time(time_count).freq(f).z_all))
                    min_z = min(min(stats(network_count).time(time_count).freq(f).z_all));
                end
                
                % Do significance-testing to get a patient-level t
                % statistic
                t_all = nan(size(ns_all,2),1);
                p_all = nan(size(ns_all,2),1);
                for tt = 2:size(ns_all,2)
                    [~,p,~,stats1] = ttest(ns_all(:,1),ns_all(:,tt));
                    t_all(tt) = stats1.tstat;
                    p_all(tt) = p;
                    
                end
                stats(network_count).time(time_count).freq(f).t_all(i,:) = t_all;
                stats(network_count).time(time_count).freq(f).p_all(i,:) = p_all;
                
            end
            
            
        end
        
    end
 
end

%% Plot
% Initialize figure
%{
nfreq (coherence) + 1 (simple) columns and 2 (2 time scales) rows
%}
figure
set(gcf,'position',[1 100 1399 500])
[ha, pos] = tight_subplot(time_count, n_freq_abs+1, [0.01 0.01], [0.1 0.05], [0.05 0.01]);

for n = 1:network_count

    net_name = stats(n).name;
    
    for t = 1:time_count
        
        % change times for x axis
        nchunks = size(stats(n).time(t).freq(1).z_all,2);
        if isfield(stats(n).time(t),'index_windows') && isempty(stats(n).time(t).index_windows) == 0
            temp_times = stats(n).time(t).index_windows(:,1)/stats(n).time(t).fs-3;
            times = realign_times(temp_times,surround_time);
        else

            times = realign_times(nchunks,surround_time);
        end
        times = round(times*1e2)/(1e2); %rounding errors
        
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
            sp = (n_freq_abs+1)*(t-1) + f + column_add;
            axes(ha(sp));
          
            z_curr = stats(n).time(t).freq(f).z_all;
            % loop over patients and plot
            for i = 1:size(z_curr,1)
                plot(times,z_curr(i,:),'ko');
                hold on
            end
            
            % plot mean F across patients
            for tt = 2:size(z_curr,2)
                
                curr_t_stats = stats(n).time(t).freq(f).t_all(:,tt);
                [~,comb_p] = ttest(curr_t_stats);
                text_out = get_asterisks(comb_p,(nchunks-1)*(n_freq_abs+1)); % should I also adjust by nfreq?
                
              %  if n == 1 && f == 5 && t == 1 && tt == 7, error('look\n'); end
                tw = stats(n).time(t).time_window;
                if strcmp(text_out,'') == 1
                    plot([times(tt)-tw/2 times(tt)+tw/2],...
                    [nanmean(z_curr(:,tt)) nanmean(z_curr(:,tt))],...
                    'k','linewidth',4);
                else
                    plot([times(tt)-tw/2 times(tt)+tw/2],...
                    [nanmean(z_curr(:,tt)) nanmean(z_curr(:,tt))],...
                    'g','linewidth',4);
                end
            end
            

            
            if t == 2 && f == 4
                 xlabel('Time relative to spike peak (s)')
            end 
           
            if t == 1 && strcmp(net_name,'coherence') == 1
                title(sprintf('%s',...
                    strrep(stats(n).time(t).freq(f).name,'_',' ')))
            elseif t == 1 && strcmp(net_name,'simple') == 1
                title('correlation')
            end
            
        end
        
        
        
    end
end

for sp = 1:length(ha)
    axes(ha(sp))
    % formatting
    ylim([min_z-1 max_z+1]);
    if mod(sp,n_freq_abs+1) == 1
        ylabel(sprintf('%s s time window',stats(1).time(floor((sp/n_freq_abs+1))).name));
    else
        yticklabels([])
    end
    
    set(gca,'fontsize',15)
end

print(gcf,[out_folder,'ns'],'-depsc');



end