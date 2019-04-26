function calc_network_metrics(whichPts)

%% Parameters
which_freq = 1; % 1 = alpha/theta; 2 = beta, 3 = low gamma, 4 = high gamma, 5 = ultra high, 6 = broadband
n_times = 11;
spike_window_times = [-0.2 0.8];

%% Get file locations, load spike times and pt structure
locations = spike_network_files;
main_folder = locations.main_folder;
results_folder = [main_folder,'results/'];
data_folder = [main_folder,'data/'];
script_folder = locations.script_folder;
addpath(genpath(script_folder));
spike_times_file = [data_folder,'spike_times/times.mat'];
pt_file = [data_folder,'spike_structures/pt.mat'];
bct_folder = locations.BCT;
addpath(genpath(bct_folder));

times = load(spike_times_file); % will result in a structure called "out"
times = times.out;
pt = load(pt_file); % will create a structure called "pt"
pt = pt.pt;

if isempty(whichPts) == 1
    whichPts = 1:length(times);
end

for whichPt = whichPts
    
    %% Prep patient
    
    % Skip it if name is empty
    if isempty(times(whichPt).name) == 1, continue; end
    name = times(whichPt).name;
    fprintf('\nDoing %s\n',name);
    pt_folder = [results_folder,name,'/'];
    adj_folder = [results_folder,name,'/adj/'];
    fs = pt(whichPt).fs;
    
    %% Load adjacency matrices and calculate metrics
    listing = dir([adj_folder,'adj*.mat']);
    if length(listing) == 0
        fprintf('No adjacency matrices for %s, skipping...\n',name);
        continue
    end
    
    % Number of spikes expected
    n_spikes = sum(~isnan(times(whichPt).spike_times));
    
    % Prep network matrices
    ge = nan(n_spikes,n_times);
    ns_seq = nan(n_spikes,n_times);
    sync = nan(n_spikes,n_times);
    ec_seq = nan(n_spikes,n_times);
    dev = nan(n_spikes,fs*n_times);
    bin_dev = nan(n_spikes,n_times);
    
    % Initialize spike count
    s_count = 0;
    
    for f = 1:length(listing)
        fname = listing(f).name;
        meta = load([adj_folder,fname]); %outputs a struct named meta
        
        % Load the adjacency matrix file
        meta = meta.meta;
        
        % Load the spike file
        spike = load([pt_folder,sprintf('spikes_%d.mat',f)]);
        spike = spike.spike;
        
        % Loop through spikes
        for s = 1:length(meta.spike)

            s_count = s_count + 1;
            
            fprintf('Doing spike %d of %d\n',s_count,length(listing)*100);
            
            adj_all_t= meta.spike(s).adj(which_freq).adj;
            seq_chs = meta.spike(s).is_seq_ch; % binary array
            sp_ch = meta.spike(s).is_sp_ch; % binary array
            
            %% Get spike signal
            data = load([pt_folder,'basic_info.mat']); % returns a structure called data
            data = data.data;
            
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
            
            % Get the values
            values = spike(s).values(:,~ignore);
            values = values(:,sp_ch);
            
            % Divide it up into windows
            peak = round(size(values,1)/2);
            index_windows = zeros(n_times,2); %6 is the spike window
            tick_window = 1*fs;
            spike_window = peak + spike_window_times*fs;
            
            for i = 1:size(index_windows,1)
                
                % if we're at the spike window (6), then the index window
                % is the spike window. If i is 1, then it's 5 seconds
                % before. If i is 11, it's 5 seconds after
                index_windows(i,:) = spike_window + tick_window*(i-6);
            end
            old_values = values;
            values = values(round(index_windows(1,1)):round(index_windows(end,2)));
            
            % readjust size
            if length(values) < size(dev,2)
                values = [values;nan(size(dev,2)-length(values),1)];
            elseif length(values) > size(dev,2)
                values = values(1:end-(length(values)-size(dev,2)));
            end
            
            dev_t = sqrt((values-nanmedian(values)).^2);
            dev(s_count,:) = dev_t;
            
            for i = 1:size(index_windows,1)
                new_dev = sqrt((old_values-nanmedian(values)).^2);
                bin_dev(s_count,i) = nanmean(new_dev(round(index_windows(i,1)):...
                    min(round(index_windows(i,2)),length(new_dev))));
            end
            
            % Loop through times
            for tt = 1:size(adj_all_t,1)
            
                % Get adj matrix of interest
                adj = squeeze(adj_all_t(tt,:,:));

                %% Calculate metrics
                ge(s_count,tt) = efficiency_wei(adj,0); % global efficiency
                sync(s_count,tt) = synchronizability_sp(adj);
                
                % node strength of all ch in seq
                ns_temp = strengths_und(adj); 
                ns_seq(s_count,tt) = sum(ns_temp(seq_chs));
                
                % eigenvector centrality of all ch in seq
                ec_temp = eigenvector_centrality_und(adj);
                ec_seq(s_count,tt) = sum(ec_temp(seq_chs));
                
            end
            
        end
        
    end
    
    %% Aggregate metrics
    avg_ns_seq = nanmean(ns_seq-median(ns_seq,2),1);
    avg_ge = nanmean(ge - median(ge,2),1);
    avg_sync = nanmean(sync - median(sync,2),1);
    avg_ec_seq = nanmean(ec_seq - median(ec_seq,2),1);
    
    %% Get avg deviation of signal
    avg_dev = nanmean(dev,1);
    avg_bin_dev = nanmean(bin_dev,1);
    
    %% Plot aggregated metrics
    figure
    subplot(2,4,1)
    plot(avg_ns_seq,'s-')
    title('Node strength of sequence chs');
    
    subplot(2,4,2)
    plot(avg_ec_seq,'s-')
    title('Eigenvector centrality of sequence chs');
    
    subplot(2,4,3)
    plot(avg_ge,'s-')
    title('Global efficiency');
    
    subplot(2,4,4)
    plot(avg_sync,'s-')
    title('Synchronizability');
    
    subplot(2,4,5)
    plot(avg_dev)
    title('Avg deviation of signal for sp channel');
    
    subplot(2,4,6)
    plot(avg_bin_dev,'s-')
    title('Avg deviation of signal for sp channel, binned');
    
   
    
    
    
    
    [p,h,stats] = signrank(sync(:,4),sync(:,7));
   
    
end


end