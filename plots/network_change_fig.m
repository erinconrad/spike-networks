function network_change_fig(which_times,remove_sig_sd,not_a_spike)

%{
I am not sure how to plot the NBS statistics, not sure what to show other
than a p-value which is kind of lame. Thinking I will just show the
permanova statistic.
%}

%% Parameters
only_simple = 0;
bonferroni_sd = 0;
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
    
    time_text = listing(i).name;
    time_window = str2num(time_text);
    time_window_folder = [sig_dev_folder,time_text,'/'];
    
    % load the file
    sub_listing = dir(time_window_folder);
    for k = 1:length(sub_listing)
        
        if contains(sub_listing(k).name,'.mat') == 0, continue; end
        
        if not_a_spike == 0
            if contains(sub_listing(k).name,'not_spike') == 1, continue; end
        else
            if contains(sub_listing(k).name,'not_spike') == 0, continue; end
        end
        count = count+1;
        temp_sig_dev = load([time_window_folder,sub_listing(k).name]);
        sig_dev(count).name = time_text;
        sig_dev(count).time_window = time_window;
        sig_dev(count).sig_dev = temp_sig_dev.sig_dev; 
    end
end

%% Now combine t statistics across patients to find the times of significant power change
%{
These t statistics are from a paired t test comparing the power in the
first time window against subsequent time windows
%}
n_windows = count;
for t = 1:n_windows
    sig_dev(t).t_stat_all = nan(length(sig_dev(1).sig_dev),...
        length(sig_dev(1).sig_dev(1).stats));
    for i = 1:length(sig_dev(t).sig_dev)
        for tt = 1:length(sig_dev(t).sig_dev(i).stats)
            if isempty(sig_dev(t).sig_dev(i).stats(tt).tstat) == 0
                % take negative so positive if later time larger
                sig_dev(t).t_stat_all(i,tt) = -sig_dev(t).sig_dev(i).stats(tt).tstat;
            end
        end
    end
    
    sig_dev(t).p_all = nan(size(sig_dev(t).t_stat_all,2),1);
    sig_dev(t).sig = zeros(size(sig_dev(t).t_stat_all,2),1);
    for tt = 2:size(sig_dev(t).t_stat_all,2)
        
        % One sample ttest on the t-statistics across patients
        [~,p] = ttest(sig_dev(t).t_stat_all(:,tt));
        sig_dev(t).p_all(tt) = p;
        if bonferroni_sd == 1
            adjusted_alpha = alpha/(size(sig_dev(t).t_stat_all,2)-1);
        else
            adjusted_alpha = alpha;
        end
        
        % If power in this time period different from first power, then say
        % this time and all subsequent times have a significant power
        % change and should be excluded
        if p < adjusted_alpha
            sig_dev(t).sig(tt:end) = 1;
            break
        end
    end
    sig_dev(t).sig(1) = 1; % Also exclude the first time
    sig_dev(t).sig = logical(sig_dev(t).sig); % indices of significant power changes
    sig_dev(t).times = round(sig_dev(t).sig_dev(1).time_window'*1e2)/(1e2); % get times
end


%% Loop through network types
listing = dir(perm_folder);
network_count = 0;
n_freq_abs = 0;
max_F = 0;

for l = 1:length(listing)
    name= listing(l).name;
    
    % Skip if . or ..
    if strcmp(name,'.') == 1 || strcmp(name,'..') == 1
        continue
    end
    
    % Skip if not a directory
    if listing(l).isdir == 0, continue; end
    
    if only_simple == 1
        if contains(name,'coherence') == 1,continue;end
    end
    
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
        
        % Skip if not one of the ones I asked for
        if ismember(time_window,which_times) == 0, continue; end
        
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
        
        
       if isfield(sim,'index_windows')
           

           
           for f = 1:nfreq
                stats(network_count).time(time_count).index_windows = sim(f).index_windows;
                stats(network_count).time(time_count).fs = sim(f).fs;
               
                stats(network_count).time(time_count).freq(f).F_all = ...
                    nan(length(pt_listing),size(sim(f).index_windows,1));

                

            end
       end
        
        
        
        pt_names = {};
        
        % loop through pts
        for i = 1:length(pt_listing)
            
            if contains(pt_listing(i).name,'not'),continue;end
            
            pt_name = pt_listing(i).name;
            pt_name_pt = strsplit(pt_name,'_');
            pt_name_pt = pt_name_pt{1};
            
            [a,b] = ismember(pt_name_pt,pt_names);
            if a == 1
                pt_idx = b;
            else
                pt_names = [pt_names;pt_name_pt];
                pt_idx = length(pt_names);
            end
            
            % load pt file
            sim = load([time_folder,pt_name]);
            sim = sim.sim;
            
            
            for f = 1:nfreq
                stats(network_count).time(time_count).freq(f).name = sim(f).name;
                if isfield(stats(network_count).time(time_count).freq(f),'F_all') && ...
                        size(stats(network_count).time(time_count).freq(f).F_all,2) > length(sim(f).F)
                    sim(f).F = [sim(f).F; ...
                        nan(size(stats(network_count).time(time_count).freq(f).F_all,2)-length(sim(f).F),1)];
                    sim(f).p = [sim(f).p; ...
                        nan(size(stats(network_count).time(time_count).freq(f).p_all,2)-length(sim(f).p),1)];
                    fprintf('\nWarning, had to pad F stats with nans for %s.\n',pt_name);
                    % this is because I originally calculated adjacency
                    % matrices for -3 to +3 for some patients, but I do
                    % later patients and subsequent analyses on just -3
                    % to 0
                end
                stats(network_count).time(time_count).freq(f).F_all(pt_idx,:) = sim(f).F;
                stats(network_count).time(time_count).freq(f).p_all(pt_idx,:) = sim(f).p;
                
                F_curr = stats(network_count).time(time_count).freq(f).F_all;
                if max_F < max(max(F_curr))
                    max_F = max(max(F_curr));
                end
                
                
                
                % convert F stats to z scores to compare across time points
                stats(network_count).time(time_count).freq(f).z_all(pt_idx,:) = (sim(f).F-nanmean(sim(f).F))./nanstd(sim(f).F);
                            
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
set(gcf,'position',[1 100 1399 600])
%[ha, pos] = tight_subplot(1, n_freq_abs+1, [0.01 0.01], [0.12 0.07], [0.05 0.01]);
if only_simple == 1
    [ha, pos] = tight_subplot(time_count, 1, [0.1 0.01], [0.12 0.07], [0.05 0.01]);
else
    [ha, pos] = tight_subplot(time_count, n_freq_abs+1, [0.1 0.01], [0.12 0.07], [0.05 0.01]);
end

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
        
        % Fix rounding error
        times = round(times*1e2)/(1e2);
        
        % Find the corresponding sig_dev time window index
        for tt = 1:length(sig_dev)
            if strcmp(sig_dev(tt).name,stats(n).time(t).name) == 1
                sig_dev_idx = tt;
                break
            end
        end

        % The times in the signal deviation structure should line up with
        % the times in the network structure
        sig_power_change_bin = sig_dev(sig_dev_idx).sig;
        sig_dev_times = sig_dev(sig_dev_idx).times;
        if ~isequal(sig_dev_times,times)
            error('Non-aligning times');
        end
        
        nfreq = length(stats(n).time(t).freq);
        
        
        if remove_sig_sd
            % Only include times without sig power change
            times = times(~sig_power_change_bin);
        end
        
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
          
            z_curr = stats(n).time(t).freq(f).z_all;
            if remove_sig_sd
                % Only include times without sig power change
                z_curr = z_curr(:,~sig_power_change_bin);
            end
            
            % loop over patients and plot
            for i = 1:size(z_curr,1)
                plot(times,z_curr(i,:),'ko');
                hold on
            end
                        
            % plot mean F across patients
            for tt = 1:size(z_curr,2)
                
                curr_p_vals = stats(n).time(t).freq(f).p_all(:,tt);
                comb_p = fisher_p_value(curr_p_vals);
                
                if remove_sig_sd == 1
                    adjusted_alpha2 = sum(~sig_power_change_bin)*(n_freq_abs+1);
                else
                    adjusted_alpha2 = (nchunks-1)*(n_freq_abs+1);
                end
                text_out = get_asterisks(comb_p,adjusted_alpha2);
                %if f == 8, error('look\n'); end
                
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
                %}
            end
            
            %if f == 7, error('look\n'); end
            
            if t == 3 && f == 4
                 xlabel('Time relative to spike peak (s)')
            end 
           
            if t == 1 && strcmp(net_name,'coherence') == 1
                title(sprintf('%s',...
                    strrep(stats(n).time(t).freq(f).name,'_',' ')))
            elseif t == 1 && strcmp(net_name,'simple') == 1
                title('correlation')
            end
            
          %  xlim([-3 0]) 
          xl = get(gca,'xlim');
          xl(2) = 0.1;
          set(gca,'xlim',xl);
            
        end
        
        
        
    end
end

for t = 1:time_count
    max_z = 0;
    min_z = 0;
    for n = 1:network_count
        nfreq = length(stats(n).time(t).freq);
        for f = 1:nfreq
            if max_z < max(max(stats(n).time(t).freq(f).z_all))
                max_z = max(max(stats(n).time(t).freq(f).z_all));
            end

            if min_z > min(min(stats(n).time(t).freq(f).z_all))
                min_z = min(min(stats(n).time(t).freq(f).z_all));
            end
        end
    end
    ylim([min_z max_z]);
end

for sp = 1:length(ha)
    axes(ha(sp))
    % formatting
    %ylim([min_z-1 max_z+1]);
    if mod(sp,n_freq_abs+1) == 1
        %ylabel(sprintf('%s s time window',stats(1).time(floor((sp/n_freq_abs+1))).name));
        if ceil(sp/(n_freq_abs+1)) == 2
            ylabel('z-score')
        end
    else
        yticklabels([])
    end
    
    if sp <= (n_freq_abs+1)*(time_count-1)
      %  xticklabels([])
    end
    
    set(gca,'fontsize',20)
end

%{
annotation('textbox',[0.01 0.78 0.2 0.2],'String','A','Fontsize',30,...
    'linestyle','none');
annotation('textbox',[0.01 0.5 0.2 0.2],'String','B','Fontsize',30,...
    'linestyle','none');
annotation('textbox',[0.01 0.19 0.17 0.2],'String','C','Fontsize',30,...
    'linestyle','none');
%}

print(gcf,[out_folder,'network_change'],'-depsc');

%% Say the patients with significant pre-spike rise
midpoint = nchunks/2;
for n = 1:network_count
    fprintf('\n for %s:\n\n',listing(n).name);
for t = 1:(time_count)
    fprintf('\n for time window %s:\n\n',time_listing(t).name);
    for i = 1:length(pt_names)
        fprintf('\n%s had significant pre-spike network change for:',pt_names{i});
        for f = 1:nfreq
            for tt = 1:midpoint - 1
            
                % Get the p-value
                p = stats(n).time(t).freq(f).p_all(i,tt);
                
                if p < 0.05/(n_freq_abs+1)/length(pt_names)/(nchunks-1)
                    fprintf('\n%s time %d.\n',freq_names{f},tt);
                end

            end
        end
        fprintf('\n\n\n');
    end
end
end

end