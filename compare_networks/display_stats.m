function display_stats(simple,time_window)

%{
simple script to display summary stats

add code to find patients for whom:
1) there is no network change at the spike (when there is a significant
signal deviation)
2) there is a pre- or post-spike network change
%}

%% Parameters
alpha = 0.05;

locations = spike_network_files;
main_folder = locations.main_folder;
results_folder = [main_folder,'results/'];
time_text = sprintf('%1.1f/',time_window);
if simple == 1
    out_folder = [results_folder,'stats_sum/simple/',time_text];
elseif simple == 0
    out_folder = [results_folder,'stats_sum/coherence/',time_text];
end

all_tables = load([out_folder,'stats_sum.mat']);
all_tables = all_tables.all_tables;

no_signal_dev_text = {};
no_spike_change_text = {};
non_spike_change_text = {};

no_signal_dev = [];
no_spike_change = [];
non_spike_change = [];

for i = 1:length(all_tables)
for j = 1:length(all_tables(i).freq)
    fprintf('%s %s\n',all_tables(i).name,all_tables(i).freq(j).name)
    all_tables(i).freq(j).table
    
    %% Check if no signal deviation
    sig_dev = all_tables(i).freq(j).table.SignalDev;
    sig_dev = cellfun(@str2num,sig_dev);
    sig_dev_not_first_time = sig_dev(2:end);
    num_sig = sum(sig_dev_not_first_time < alpha/length(sig_dev_not_first_time)); % get number of significant signal deviations
    if num_sig == 0
        % If no significant signal deviations, add the name and frequency
        % to table
        no_signal_dev_text = [no_signal_dev_text;[all_tables(i).name,' ',all_tables(i).freq(j).name]];
        no_signal_dev = [no_signal_dev; i,j];
    else
        significant_sig_dev = find(sig_dev_not_first_time < alpha/length(sig_dev_not_first_time));
    end
        
    %% Now check for lack of network deviation at times of signal deviation
    any_sig_perm_or_nbs = 0;
    perm = cell2mat(cellfun(@str_2_num_asterisk,all_tables(i).freq(j).table.Perm,'UniformOutput',false));
    nbs = cell2mat(cellfun(@str_2_num_asterisk,all_tables(i).freq(j).table.NBS,'UniformOutput',false));
    
    for k = 1:length(significant_sig_dev)
        % add 1 to the row (because it's an array excluding the first time)
        sig_time = significant_sig_dev(k) + 1;
        
        % get the corresponding permutation and nbs p values
        perm_time = perm(sig_time);
        nbs_time = nbs(sig_time);
        
        if perm_time < alpha/(length(perm)-1) || nbs_time < alpha/(length(perm)-1)
            any_sig_perm_or_nbs = 1;
            break
        end
        
    end
    % If never found a significant perm or nbs at any of the signal
    % deviation times
    if any_sig_perm_or_nbs == 0
        no_spike_change_text = [no_spike_change_text;[all_tables(i).name,' ',all_tables(i).freq(j).name]];
        no_spike_change = [no_spike_change;i,j];
    end
    
    %% Now check for any network deviations not at times of signal deviation
    significant_sig_dev_times = significant_sig_dev + 1; % add 1 
    all_times = 1:length(all_tables(i).freq(j).table.NBS);
    no_sig_dev_times = all_times(~ismember(all_times,significant_sig_dev_times));
    
  %  if strcmp(all_tables(i).name,'HUP078') == 1, error('look\n'); end
    
    % Loop through these times without significant signal dev
    for k = 1:length(no_sig_dev_times)
        curr_time = no_sig_dev_times(k);
        perm_time = perm(curr_time);
        nbs_time = nbs(curr_time);
        
        if perm_time < alpha/(length(perm)-1)
            non_spike_change_text = [non_spike_change_text;[all_tables(i).name,' ',...
                all_tables(i).freq(j).name,' ',curr_time,' perm']];
            non_spike_change = [non_spike_change;i,j,(curr_time),1];
        end
        
        if nbs_time < alpha/(length(perm)-1)
            non_spike_change_text = [non_spike_change_text;[all_tables(i).name,' ',...
                all_tables(i).freq(j).name,' ',curr_time,' nbs']];
            non_spike_change = [non_spike_change;i,j,(curr_time),2];
        end
    end
    
    %pause
end
end

%% Summarize results
method = {'perm','nbs'};
fprintf(['\nThe following are patients for whom we did not detect a signficant\n'...
    'change in signal deviation at any time:\n']);
for i = 1:size(no_signal_dev,1)
    fprintf('%s\n',all_tables(no_signal_dev(i,1)).name);
end

fprintf(['\nThe following are patients for whom we did not detect a signficant\n'...
    'change in network at any of the signal deviation times:\n']);
n_freq = length(all_tables(1).freq);
% find the patients for whom there are n_freq occurrences (implying no
% change at any frequency)
unique_pts = unique(no_spike_change(:,1));
for i = 1:length(unique_pts)
    flag = 0;
    curr_pt = unique_pts(i);
    n_pt_freq = sum(curr_pt == no_spike_change(:,1)); % number of frequencies for that patient
    if n_freq == n_pt_freq
        % if the number of frequencies without significant network change
        % for that patient is the total number of frequencies, flag the
        % patient
        flag = 1;
    end
    if flag == 1
        fprintf('%s\n',all_tables(curr_pt).name);
    end
end

fprintf(['\nThe following are patients, frequencies, times, and method for whom we detected\n'...
    'a significant network change during a time of no signal deviation:\n']);
for i = 1:size(non_spike_change,1)
    fprintf('Pt: %s, freq: %s, time: %1.1f, method %s\n',...
        all_tables(non_spike_change(i,1)).name,...
        all_tables(non_spike_change(i,1)).freq(non_spike_change(i,2)).name,...
        all_tables(non_spike_change(i,1)).freq(non_spike_change(i,2)).table.Time(non_spike_change(i,3)),...
        method{non_spike_change(i,4)});
end

end

function y = str_2_num_asterisk(x)
C = strsplit(x,'*');
y = C{1};
y = str2double(y);
end
