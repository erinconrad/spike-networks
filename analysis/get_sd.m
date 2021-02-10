function sig_dev = get_sd


do_avg_text = '';

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
%sig_dev_folder = [results_folder,'signal_deviation/manual/',do_avg_text];
sig_dev_folder = [results_folder,'power/manual/'];
perm_folder = [results_folder,'perm_stats/'];
ers_folder = [results_folder,'ers/'];
ns_folder = [results_folder,'metrics/manual/'];

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
    
    if contains(listing(i).name,'avg_chs'), continue; end
    
    time_text = listing(i).name;
    time_window = str2num(time_text);
    time_window_folder = [sig_dev_folder,time_text,'/'];
    
    count = count+1;
    
    % load the file
    sub_listing = dir(time_window_folder);
    for k = 1:length(sub_listing)
       
        if contains(sub_listing(k).name,'.mat') == 0, continue; end
        

        temp_sig_dev = load([time_window_folder,sub_listing(k).name]);

        
        if contains(sub_listing(k).name,'not') == 1
            sig_dev(count).is_spike(2).sig_dev = temp_sig_dev.sig_dev; 
            sig_dev(count).is_spike(2).name = 'not_spike';
            sig_dev(count).is_spike(2).name = time_text;
            sig_dev(count).is_spike(2).time_window = time_window;
        else
            sig_dev(count).is_spike(1).sig_dev = temp_sig_dev.sig_dev; 
            sig_dev(count).is_spike(1).name = 'spike';
            sig_dev(count).is_spike(1).name = time_text;
            sig_dev(count).is_spike(1).time_window = time_window;
        end
    end
end

n_windows = count;
for t = 1:n_windows % loop over time windows
    for s = 1:2 % loop over spike and not a spike
        sig_dev(t).is_spike(s).times = round(sig_dev(t).is_spike(s).sig_dev(1).time_window'*1e2)/(1e2); % get times
    end
end


%% Now combine t statistics across patients to find the times of significant power change
%{
These t statistics are from a paired t test comparing the power in the
first time window against subsequent time windows
%}
%{
n_windows = count;
for t = 1:n_windows % loop over time windows
    for s = 1:2 % loop over spike and not a spike
        
        if s>length(sig_dev(t).is_spike), continue; end
        sig_dev(t).is_spike(s).t_stat_all = nan(length(sig_dev(t).is_spike(s).sig_dev),...
            length(sig_dev(t).is_spike(s).sig_dev(1).stats));
        for i = 1:length(sig_dev(t).is_spike(s).sig_dev)
            for tt = 1:length(sig_dev(t).is_spike(s).sig_dev(i).stats)
                if isempty(sig_dev(t).is_spike(s).sig_dev(i).stats(tt).tstat) == 0
                    % take negative so positive if later time larger
                    sig_dev(t).is_spike(s).t_stat_all(i,tt) = -sig_dev(t).is_spike(s).sig_dev(i).stats(tt).tstat;
                end
            end
        end

        sig_dev(t).is_spike(s).p_all = nan(size(sig_dev(t).is_spike(s).t_stat_all,2),1);
        sig_dev(t).is_spike(s).sig = zeros(size(sig_dev(t).is_spike(s).t_stat_all,2),1);

        sig_dev(t).is_spike(s).times = round(sig_dev(t).is_spike(s).sig_dev(1).time_window'*1e2)/(1e2); % get times

    end
end

%% Decide what is a significant power change
for t = 1:n_windows % loop over time windows
    
    if bf == 1
        adjusted_alpha = alpha/(size(sig_dev(t).is_spike(1).t_stat_all,2)-1);
    else
        adjusted_alpha = alpha;
    end
    sig_dev(t).tests.unpaired.sig = zeros(size(sig_dev(t).is_spike(1).t_stat_all,2),1);
    sig_dev(t).tests.paired.sig = zeros(size(sig_dev(t).is_spike(1).t_stat_all,2),1);
    for tt = 2:size(sig_dev(t).is_spike(1).t_stat_all,2)
    
        % Do unpaired test of tstats
        [~,p] = ttest(sig_dev(t).is_spike(1).t_stat_all(:,tt));
        if p < adjusted_alpha
            sig_dev(t).tests.unpaired.sig(tt:end) = 1;
            
        end
        
        if s > length(sig_dev(t).is_spike), continue; end
        if length(sig_dev(t).is_spike(1).t_stat_all(:,tt)) ~=...
            length(sig_dev(t).is_spike(2).t_stat_all(:,tt))
            continue;
        end

        % Do paired test of tstats, comparing against not a spike
        [~,p] = ttest(sig_dev(t).is_spike(1).t_stat_all(:,tt),...
            sig_dev(t).is_spike(2).t_stat_all(:,tt));
        if p < adjusted_alpha
            sig_dev(t).tests.paired.sig(tt:end) = 1;
            
        end
    
    end
    sig_dev(t).tests.unpaired.sig(1) = 1; % also exclude first time
    sig_dev(t).tests.paired.sig(1) = 1;
    
end
%}


end