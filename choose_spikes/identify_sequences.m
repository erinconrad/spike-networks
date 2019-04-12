function out = identify_sequences(pt,cluster)

%{
- this function will take the pt and cluster structures generated as part
of the spike project and will randomly select 1,000 spikes/not spikes for 
each patient and will record the times
%}

%% Parameters
n_spikes = 1e3; % how many spikes to select
n_skip = 1e3; % skip patient if <1000 sequences
window_sp = 11; % don't allow sequences closer than 11 seconds apart
window_ns = 5; % Don't allow not-spikes within 6 seconds of spikes
range_time_to_add = 60; % will look for not_spikes within 6 to 60 seconds after spikes
removeTies = 1;
time_attempt = 120; % allow < 1000 spikes if takes more than 120 seconds

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
    
    %% Patient parameters
    szTimes = pt(whichPt).newSzTimes;    
    seq_matrix = pt(whichPt).seq_matrix;
    
    %% Remove ties
    if removeTies == 1
        keep = ones(size(seq_matrix,2),1);
        for s = 1:size(seq_matrix,2)
           curr_seq = seq_matrix(:,s);
           nonans = curr_seq(~isnan(curr_seq));
           norepeats = unique(nonans);
           if length(norepeats) < 0.5*length(nonans)
               keep(s) = 0;
           end
        end
        seq_matrix(:,keep==0) = [];
        fprintf(['%s had %d sequences (%1.2f of all sequences) deleted'...
        'for having >50 percent ties\n%d sequences remain\n'],...
        pt(whichPt).name,sum(keep == 0),sum(keep == 0)/length(keep),sum(keep==1));

    end
    
    %% Remove ictal sequences
    all_times = seq_matrix(:);
    icTimes = find(any(all_times >= (szTimes(:,1)-repmat(60,size(szTimes,1),1))' ...
        & all_times <= szTimes(:,2)',2));
    seq_matrix(icTimes) = nan;
    fprintf('Removed %d ictal spikes\n',length(icTimes));
    
    %% Get cluster info
    all_times_all = cluster(whichPt).all_times_all; % all spike times
    all_spikes = cluster(whichPt).all_spikes; % all spike channels
    idx = cluster(whichPt).idx; % the cluster index for every spike
    bad_cluster = cluster(whichPt).bad_cluster; % which clusters are bad
    
    %% Compare number of spikes in cluster array and my data
    if sum(sum(~isnan(seq_matrix))) ~= length(all_times_all)
        error('Warning, number of spikes do not align\n');
    end
    
 
    %% Find bad spikes
    bad_idx = find(ismember(idx,bad_cluster));
    
    % Nx2 array of bad spikes, showing the channel and time
    bad_spikes = [all_spikes(bad_idx),all_times_all(bad_idx)];
    
    %% Get all sequences   
    new_seq_matrix = seq_matrix;
    n_removed = 0;
    
    %% Go through sequence matrix and remove bad spikes
    for ich = 1:size(seq_matrix,1)
        % loop across electrodes
        
        % All spike times for this channel
        spikeTimesCh = seq_matrix(ich,:);
        
        % Get the bad spikes in that channel
        bad_times_for_ch = bad_spikes(bad_spikes(:,1) == ich,2);
        
        % Make sure I am finding all of them
        Lia = ismember(spikeTimesCh,bad_times_for_ch);
        if sum(Lia) ~= length(bad_times_for_ch)
            error(sprintf('Did not find all bad spikes for channel %d\n',ich));
        end

        n_removed = n_removed + sum(Lia);
        
        % Make bad spikes nans
        spikeTimesCh(Lia==1) = nan;
        new_seq_matrix(ich,:) = spikeTimesCh;
        
        
    end
    
    if n_removed~=length(bad_idx)
        error('Incorrect number of bad spikes removed\n');
    end
    
    seq_matrix = new_seq_matrix;
    fprintf('%d sequences remain\n',size(seq_matrix,2));
    
    % I now have a bunch of things that I'm pretty sure are spikes
    %% Get first channel and time of each sequence
    [first_time,first_ch] = min(seq_matrix,[],1);
    
    % Skip if < 1,000 sequences
    if length(first_time) < n_skip
        fprintf('Skipping %s due to not enough sequences.\n',name);
        continue
    else
        fprintf('%s has %d sequences, proceeding\n',name,length(first_time));
    end
    
    % Check if enough sequences far apart
    all_time_diff = abs(diff(first_time));
    n_f_apart = find(all_time_diff > window_sp);
    if length(n_f_apart) < n_spikes
        fprintf('Skipping %s due to not enough far apart sequences.\n',name);
        continue
    else
        fprintf('%s has %d sequences far enough apart, proceeding\n',name,length(n_f_apart));
    end
    
    %% Get spike times
    tic
    for i = 1:length(spike_times)
        
        while 1
            
            % Pick random spike time
            ix = randi(length(first_time));
            t = first_time(ix);
            ch = first_ch(ix);
            
            % don't allow it if it's an ignored channel per the json file
            label = pt(whichPt).electrodeData.electrodes(ch).name;
            if ismember(label,pt(whichPt).ignore.names) == 1
                continue;
            end
            
            % Loop through prior spike times and see if it is too close
            too_close = 0;
            for j = 1:i-1
                
                if abs(t-spike_times(j)) < window_sp
                    too_close = 1;
                end
                
            end
            
            % If it is not too close to any prior spike, now look to see 
            % if it is too close to other spike sequences
            if too_close == 0
                
                time_diff_all = abs(t - first_time);
                time_diff_all = sort(time_diff_all);
                in_window = find(time_diff_all < window_sp);
                
                % If there is only one spike (itself) within allowable
                % window, add it
                if length(in_window) == 1

                    spike_times(i) = t;
                    spike_chs(i) = ch;
                    break
                end
            end
            a = toc;
            if a > time_attempt 
                break
            end
        end
        
    end
    
    %% Get not-a-spike times
    for i = 1:length(not_spike_times)
        
        while 1
            
           
            % Pick random spike
            ix = randi(length(first_time));
            
            % add random time to look after it. I am doing this because
            % perhaps there is something funny that happens N seconds after
            % every spike. By looking a random amount of time after the
            % spike, then I should be more randomly sampling non-spike
            % times
            t = first_time(ix) + randi([window_ns,range_time_to_add]);
            
            % subtract this time from all spikes times to see how close it
            % is
            time_diff = sort(abs(t - first_time));
            
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
        if isnan(spike_chs(i)) == 1, continue; end
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
    
    
    %% Show some summary plots for patient
    figure
    subplot(1,2,1)
    histogram(out(whichPt).spike_times)
    title('Spike times')

    subplot(1,2,2)
    histogram(out(whichPt).not_spike_times)
    title('Not spike times')

    pause
    close(gcf)
    
end

%% See how many we have
np = 0;
for i = 1:length(out)
    if isempty(out(i).name) == 0
        np = np+1;
    end
end

fprintf('We have %d patients\n',np);


    

end