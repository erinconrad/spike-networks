function sig_deviation(whichPts)

%{
The goal of this function is to compare time periods around the spike to
determine if the signal deviation is significantly different across time
periods. If it IS significantly different, I will ignore network
differences between those time periods.

%}

%% Parameters
nboot = 1e3;
time_window = 0.5; %in seconds
n_chunks = 22; % number of time windows


%% Get file locations, load spike times and pt structure
locations = spike_network_files;
main_folder = locations.main_folder;
results_folder = [main_folder,'results/'];
data_folder = [main_folder,'data/'];
script_folder = locations.script_folder;
addpath(genpath(script_folder));
spike_times_file = [data_folder,'spike_times/times.mat'];
pt_file = [data_folder,'spike_structures/pt.mat'];

times = load(spike_times_file); % will result in a structure called "out"
times = times.out;
pt = load(pt_file); % will create a structure called "pt"
pt = pt.pt;

if isempty(whichPts) == 1
    whichPts = 1:length(times);
end



% Loop through patients
for whichPt = whichPts
    
    %% Prep patient
    
    % Skip it if name is empty
    if isempty(times(whichPt).name) == 1, continue; end
    name = times(whichPt).name;
    fprintf('\nDoing %s\n',name);
    
     % output folder
    out_folder = [results_folder,'sig_dev/',name,'/'];

    if exist(out_folder,'dir') == 0
        mkdir(out_folder);
    end
    
    
    % Get basic data about patient
    pt_folder = [results_folder,name,'/'];
    data = load([pt_folder,'basic_info.mat']); % returns a structure called data
    data = data.data;
    
    % convert ch labels to nice form and decide which to ignore
    ch_labels = data.chLabels(:,1);
    ignore = zeros(length(ch_labels),1);
    
    for i = 1:length(ch_labels)
        ch_labels{i} = ieeg_ch_parser(ch_labels{i});
        for j = 1:length(pt(whichPt).ignore.names)
            if strcmp(pt(whichPt).ignore.names(j),ch_labels{i}) == 1
                ignore(i) = 1;
            end
        end
    end
    
    
    
    % Loop through spike files
    s_count = 0;
    
    listing = dir([pt_folder,'spikes_*.mat']);
    for f = 1:length(listing)
        
        
        spike = load([pt_folder,listing(f).name]);
        spike = spike.spike;
 
        % Loop through spikes
        for s = 1:length(spike)
            %fprintf('Doing file %d of %d...\n',f,10);
            %fprintf('Doing spike %d of %d...\n',s,length(spike));
            %tic
            if isempty(spike(s).time) == 1, continue; end
            
            s_count = s_count + 1;
            
            % Grab the appropriate channels
            values = spike(s).values(:,~ignore);
            is_sp_ch = strcmp(ch_labels(~ignore),spike(s).label);
            
            %% Figure out times for which I will be calculating adjacencies
            % The peak should be the very center of each file
            peak = round(size(values,1)/2);
            
            % Get index windows
            index_windows = zeros(n_chunks,2);
            dev_windows = zeros(n_chunks,2);
            tick_window = time_window*data.fs;
            
            for i = 1:n_chunks
                index_windows(i,1) = peak - tick_window*n_chunks/2 + tick_window*(i-1);
                index_windows(i,2) = peak - tick_window*n_chunks/2 + tick_window*(i);
                
                dev_windows(i,1) = (i-1)*tick_window+1;
                dev_windows(i,2) = i*tick_window;
            end
            
            % For first spike, initialize array
            if f==1 && s == 1
                all_dev = zeros(1000,floor(index_windows(end,2) - index_windows(1,1) +1));
            end
            
            % Get values from spike channel
            curr_values = values(index_windows(1,1):index_windows(end,2),is_sp_ch);
            
            % Get deviation from median for each point
            curr_dev = sqrt((curr_values - median(curr_values)).^2);
            
            % Fill array with deviation
            all_dev(s_count,:) = curr_dev;
            
        end
        
    end
    
    % Now get the average deviation in each time period
    bin_dev = zeros(size(all_dev,1),size(index_windows,1));
    for i = 1:n_chunks
        bin_dev(:,i) = mean(all_dev(:,dev_windows(i,1):dev_windows(i,2)),2);
    end
    
    % plot binned dev
    if 0
        figure
        plot(nanmean(bin_dev,1))
    end
    
    sig_dev.bin_dev = bin_dev;
    
    
    first_sec = bin_dev(:,1);
    for i = 2:n_chunks
        
        temp_sec = bin_dev(:,i);
        
        % Do a Wilcoxon signed rank test comparing the first second to each
        % subsequent second
        [p,~,stats] = signrank(first_sec,temp_sec);
        
        sig_dev.p(i) = p;
        sig_dev.stats(i) = stats;
        
    end
    
    save([out_folder,'sig_dev.mat'],'sig_dev')
    
    %{
    % Now do a permutation test where I randomly permute the time period
    % identities for each spike
    for ib = 1:nboot
        
        % Clever way suggested by Matlab user to generate 1,000
        % (size(all_dev,1)) rows of random permutations of
        % 1:size(index_windows,1) (1:22)
        [~, out] = sort(rand(size(all_dev,1),size(index_windows,1)),2);
        
      %  temp_bin_dev = 
        
    end
    %}
    
end

end