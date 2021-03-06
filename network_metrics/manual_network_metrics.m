function manual_network_metrics(overwrite,simple,time_window)

%% Get file locations, load spike times and pt structure
locations = spike_network_files;
main_folder = locations.main_folder;
results_folder = [main_folder,'results/'];
data_folder = [main_folder,'data/'];
script_folder = locations.script_folder;
addpath(genpath(script_folder));
addpath(genpath(locations.BCT));
pt_file = [data_folder,'spike_structures/pt.mat'];
bct_folder = locations.BCT;
addpath(genpath(bct_folder));
time_text = sprintf('%1.1f/',time_window);


% EEG data folder
eeg_folder = [results_folder,'eeg_data/'];

% Adj mat folder
if simple == 1
    adj_folder = [results_folder,'adj_mat/manual/adj_simple/',time_text];
    metrics_folder = [results_folder,'metrics/manual/simple/',time_text];
else
    adj_folder = [results_folder,'adj_mat/manual/adj_coherence/',time_text];
    metrics_folder = [results_folder,'metrics/manual/coherence/',time_text];
end

if exist(metrics_folder,'dir') == 0
    mkdir(metrics_folder);
end

listing = dir([adj_folder,'*_adj.mat']);

for i = 1:length(listing)
    
    filename = listing(i).name;
    name_sp = split(filename,'_');
    name = name_sp{1};
    
    if overwrite == 0
        if exist([metrics_folder,name,'_network_stats.mat'],'file') ~= 0
            fprintf('Already did %s, skipping...\n',name);
            continue;
        end
    end
    
    metrics.name = name;
    
    % load adj matrix
    meta = load([adj_folder,filename]);
    meta = meta.meta;
    
    % load eeg data
    spike = load([eeg_folder,name,'_eeg.mat']);
    spike = spike.spike;
    
    % get sizes for matrices
    n_chs = size(meta.spike(1).adj(1).adj,3);
    n_f = length(meta.spike(1).adj);
    n_times = size(meta.spike(1).index_windows,1);
    n_spikes = length(meta.spike);
    
    % Initialize
    ge = nan(n_f,n_spikes,n_times);
    ns = nan(n_f,n_spikes,n_times,n_chs);
    bc = nan(n_f,n_spikes,n_times,n_chs);
    
    % Initialize spike count
    s_count = 0;
    
    for s = 1:length(meta.spike)

        s_count = s_count + 1;

        fprintf('Doing spike %d of %d\n',s_count,n_spikes);
        involved = spike(s).involved;
        
        for f = 1:n_f
            adj_all_t= meta.spike(s).adj(f).adj;
            for tt = 1:size(adj_all_t,1)
                
                % Get adj matrix of interest
                adj = squeeze(adj_all_t(tt,:,:));
                
                %% Calculate metrics
                
                % global efficiency of full matrix
                ge(f,s,tt) = efficiency_wei(adj,0); 
                
                % node strength of all chs
                ns(f,s,tt,:) = strengths_und(adj); 
                
                % betweenness centrality of all chs
                bc(f,s,tt,:) = betweenness_wei(adj);
                
                
                
            end
            
        end
           
    end
    
    % Average over all spikes
    avg_ge = (nanmean(ge,2));
    avg_ns = (nanmean(ns,2));
    avg_bc = (nanmean(bc,2));
        
    % Fill structure
    for f = 1:n_f
        if isfield(meta.spike(1).adj(f),'name') == 1
            metrics.freq(f).name = meta.spike(1).adj(f).name;
        else
            metrics.freq(f).name = 'correlation';
        end
        metrics.freq(f).ge.name = 'global efficiency';
        metrics.freq(f).ge.data = squeeze(avg_ge(f,:,:));
        metrics.freq(f).ns.name = 'node strength';
        metrics.freq(f).ns.data = squeeze(avg_ns(f,:,:,:));
        metrics.freq(f).bc.name = 'betweenness centrality';
        metrics.freq(f).bc.data = squeeze(avg_bc(f,:,:,:));
    end
    
    metrics.involved = involved;
    metrics.index_windows = meta.spike(1).index_windows;
    
    save([metrics_folder,name,'_network_stats.mat'],'metrics');
    
end

end