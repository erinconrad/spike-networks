function out = auto_spike_times(pt,cluster)

%{
- this function will take the pt and cluster structures generated as part
of the spike project and will randomly select 1,000 spikes/not spikes for 
each patient and will record the times
%}

%% Parameters
n_spikes = 1e3; % how many spikes to select

%% Define which patients (all with >= 1 good cluster)
whichPts = [];
for i = 1:length(pt)
    if isempty(pt(i).seq_matrix) == 0
        if size(cluster(i).bad_cluster) < cluster(i).k
            whichPts = [whichPts,i];
        end
    end
end

for whichPt = whichPts
    
    %% Initialize arrays
    spike_times = nan(n_spikes,1);
    spike_chs = nan(n_spikes,1);
    spike_labels = cell(n_spikes,1);
    not_spike_times = nan(n_spikes,1);
    
    name = pt(whichPt).name;
    fprintf('Doing %s\n',name);
    
    %% Get patient parameters
    szTimes = pt(whichPt).newSzTimes;
    
    % Reorder seizure times if out of order
    oldSzTimes = szTimes;
    szTimes = sort(szTimes,1);
    if isequal(oldSzTimes,szTimes) == 0
        error('WARNING!!! %s seizure times out of order\n',pt(whichPt).name);
    end
    
    % Combine nearly equal seizure times
    newIdx = 2;
    newSzTimes = [];
    newSzTimes(1,:) = szTimes(1,:);
    for j = 2:size(szTimes,1)
        if abs(szTimes(j,1)-szTimes(j-1,1)) < 10 && ...
                abs(szTimes(j,2)-szTimes(j-1,2))
           newIdx = newIdx - 1; 
           newSzTimes(newIdx,1) = min(szTimes(j,1),szTimes(j-1,1));
           newSzTimes(newIdx,2) = max(szTimes(j,2),szTimes(j-1,2));  
        else
           newSzTimes(newIdx,:) = szTimes(j,:);
        end
        newIdx = newIdx + 1;
    end
    
    if isequal(newSzTimes,szTimes) == 0
        error('WARNING!!! %s had duplicate seizure times\n',pt(whichPt).name);
    end
    
    %% Get cluster info and remove bad clusters   
    all_spike_times = cluster(whichPt).all_times_all; % all spike times
    idx = cluster(whichPt).idx; % the cluster index for every spike
    all_spike_locs = cluster(whichPt).all_spikes; % spike locations
    bad_cluster = cluster(whichPt).bad_cluster; % which clusters are bad
    
    % Confirm that I do not have any ictal spikes
    t = find(any(all_spike_times >= szTimes(:,1)' & all_spike_times <= szTimes(:,2)',2));
    if isempty(t) == 0
        error('WARNING: Remaining ictal spikes for %s!\n',pt(whichPt).name);
        all_spike_times(t) = [];
        idx(t) = [];
    end
    

    % Remove bad clusters
    bad_idx = (ismember(idx,bad_cluster));
    all_spike_times(bad_idx) = [];
    n_true_spikes = length(all_spike_times);
    all_spike_locs(bad_idx) = [];
    
    % Order spike times
    all_spike_times_old = all_spike_times;
    [all_spike_times,I] = sort(all_spike_times);
    all_spike_locs = all_spike_locs(I);
    if isequal(all_spike_times,all_spike_times_old) == 0
        error('Spike times not sorted\n');
    end
    
    % I now have a bunch of things that I'm pretty sure are spikes
    
    % Skip if < 10,000 spikes
    if length(all_spike_times) < n_skip
        fprintf('Skipping %s due to not enough spikes.\n',name);
        continue
    else
        fprintf('%s has %d spikes, proceeding\n',name,length(all_spike_times));
    end
    
    %% Get spike times
    for i = 1:length(spike_times)
        
        while 1
            
            % Pick random spike time
            ix = randi(n_true_spikes);
            t = all_spike_times(ix);
            ch = all_spike_locs(ix);
            
            % Loop through prior spike times and see if it is too close
            too_close = 0;
            for j = 1:i-1
                
                if abs(t-spike_times(j)) < window_sp
                    too_close = 1;
                end
                
            end
            
            % If it is not too close to any prior spike, add this spike and
            % break the loop, moving onto the next spike. 
            if too_close == 0
                spike_times(i) = t;
                spike_chs(i) = ch;
                break
            end
            
        end
        
    end
    
    %% Get not-a-spike times
    for i = 1:length(not_spike_times)
        
        while 1
            
           
            % Pick random spike
            ix = randi(n_true_spikes);
            
            % add random time to look after it. I am doing this because
            % perhaps there is something funny that happens N seconds after
            % every spike. By looking a random amount of time after the
            % spike, then I should be more randomly sampling non-spike
            % times
            t = all_spike_times(ix) + randi([window_ns,range_time_to_add]);
            
            % subtract this time from all spikes times to see how close it
            % is
            time_diff = sort(abs(t - all_spike_times));
            
            % find time differences less than allowed amount
            too_close = time_diff < window_ns;
            
            % If none are too close, now make sure it isn't too close to
            % prior not-spikes
            if sum(too_close) == 0
                
                too_close_not_spikes = 0;
                
                for j = 1:i
                   
                    if abs(t-not_spike_times(j)) < window_sp
                        too_close_not_spikes = 1;
                    end
                    
                end
                
                % If it is also not too close to prior not-spikes, add it
                % and move on to the next not-spike
                if too_close_not_spikes == 0
                    not_spike_times(i) = t;
                    break
                end
                
            end
            
            
        end
        
    end
    
    
    %% Add spike and not spike_times
    out(whichPt).name = name;
    [out(whichPt).spike_times,I] = sort(spike_times);
    out(whichPt).not_spike_times = sort(not_spike_times);
    spike_chs = spike_chs(I);
    
    % Convert spike chs to labels
    for i = 1:length(spike_chs)
        label_temp = pt(whichPt).electrodeData.electrodes(spike_chs(i)).name;
        spike_labels{i} = label_temp;
        found_it = 0;
        for j = 1:length(pt(whichPt).new_elecs.names)
            if strcmp(pt(whichPt).new_elecs.names{j},label_temp) == 1
                found_it = 1;
            end
        end
        if found_it == 0
            fprintf('Could not find %s in new labels for %s\n',label_temp,...
                pt(whichPt).name);
        end
    end
    out(whichPt).spike_labels = spike_labels;
    
end

%% See how many we have
np = 0;
for i = 1:length(out)
    if isempty(out(i).name) == 0
        np = np+1;
    end
end

fprintf('We have %d patients\n',np);

%% Show some summary plots for patient 4
figure
subplot(2,2,1)
histogram(out(4).spike_times)
title('Spike times')

subplot(2,2,2)
histogram(out(4).not_spike_times)
title('Not spike times')

subplot(2,2,3)
histogram(diff(out(4).spike_times))
title('Time between spikes')

subplot(2,2,4)
histogram(diff(out(4).not_spike_times))
title('Time between not spikes');
    

end