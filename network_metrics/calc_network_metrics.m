function calc_network_metrics(whichPts)

%% Parameters
which_freq = 2; % 1 = alpha/theta; 2 = beta, 3 = low gamma, 4 = high gamma, 5 = ultra high, 6 = broadband
n_times = 11;

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
    adj_folder = [results_folder,name,'/adj/'];
    
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
    ns_sp = nan(n_spikes,n_times);
    ns_seq = nan(n_spikes,n_times);
    sync = nan(n_spikes,n_times);
    
    % Initialize spike count
    s_count = 0;
    
    for f = 1:length(listing)
        fname = listing(f).name;
        meta = load([adj_folder,fname]); %outputs a struct named meta
        
        % Load the adjacency matrix file
        meta = meta.meta;
        
        % Loop through spikes
        for s = 1:length(meta.spike)
            
            fprintf('Doing spike %d\n',s);
            
            s_count = s_count + 1;
            
            adj_all_t= meta.spike(s).adj(which_freq).adj;
            seq_chs = meta.spike(s).is_seq_ch; % binary array
            sp_ch = meta.spike(s).is_sp_ch; % binary array
            
            % Loop through times
            for tt = 1:size(adj_all_t,1)
            
                % Get adj matrix of interest
                adj = squeeze(adj_all_t(tt,:,:));

                %% Calculate metrics
                ge(s_count,tt) = efficiency_wei(adj,0); % global efficiency
                sync(s_count,tt) = synchronizability_sp(adj);
                
                % node strength of sp channel and all ch in seq
                ns_temp = strengths_und(adj); 
                ns_sp(s_count,tt) = ns_temp(sp_ch);
                ns_seq(s_count,tt) = sum(ns_temp(seq_chs));
                
                
            end
            
        end
        
    end
    
    %% Aggregate metrics
    avg_ns_sp = nanmean(ns_sp-median(ns_sp,2),1);
    avg_ns_seq = nanmean(ns_seq-median(ns_seq,2),1);
    avg_ge = nanmean(ge - median(ge,2),1);
    avg_sync = nanmean(sync - median(sync,2),1);
    
    
    [p,h,stats] = signrank(sync(:,1),sync(:,5));
   
    
end


end