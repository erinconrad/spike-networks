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
nbs_folder = [results_folder,'nbs_stats/'];


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
n_freq_abs = 0;
max_F = 0;
max_z = 0;
min_z = 0;
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
        if n_freq_abs < nfreq
            n_freq_abs = nfreq;
        end
        
       
        
        for f = 1:nfreq
            stats(network_count).time(time_count).freq(f).F_all = ...
                nan(length(pt_listing),surround_time*2/time_window);

            stats(network_count).time(time_count).freq(f).p_all = ...
                nan(length(pt_listing),surround_time*2/time_window);
            
            stats(network_count).time(time_count).freq(f).z_all = ...
                nan(length(pt_listing),surround_time*2/time_window);
            
        end
        
        % loop through pts
        for i = 1:length(pt_listing)
            
            pt_name = pt_listing(i).name;
            pt_name_pt = strsplit(pt_name,'_');
            pt_name_pt = pt_name_pt{1};
            
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
        nchunks = size(stats(n).time(t).freq(1).F_all,2);
        times = realign_times(nchunks,surround_time);
        
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
            for tt = 1:size(z_curr,2)
                
                curr_p_vals = stats(n).time(t).freq(f).p_all(:,tt);
                comb_p = fisher_p_value(curr_p_vals);
                text_out = get_asterisks(comb_p,(nchunks-1)*(n_freq_abs+1)); % should I also adjust by nfreq?
                
                if strcmp(text_out,'') == 1
                    plot([times(tt)-0.25 times(tt)+0.25],...
                    [mean(z_curr(:,tt)) mean(z_curr(:,tt))],...
                    'k','linewidth',4);
                else
                    plot([times(tt)-0.25 times(tt)+0.25],...
                    [mean(z_curr(:,tt)) mean(z_curr(:,tt))],...
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

print(gcf,[out_folder,'network_change'],'-depsc');

end