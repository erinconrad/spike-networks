function calculate_ns(time_window,not_spike)

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
biggest_dev_folder = [results_folder,'biggest_dev/'];
seq_folder = [results_folder,'seq_data/'];


% Adj mat folder
adj_folder = [results_folder,'adj_mat/manual/adj_coherence/',time_text];
out_folder = [results_folder,'ns/',time_text];

if exist(out_folder,'dir') == 0
    mkdir(out_folder);
end

if not_spike == 1
    not_spike_text = '_not_spike';
else
    not_spike_text = '';
end

listing = dir([adj_folder,'*_adj.mat']);

for i = 1:length(listing)
    
    filename = listing(i).name;
    name_sp = split(filename,'_');
    name = name_sp{1};
    if not_spike
        if ~contains(filename,'not'), continue; end
    else
        if contains(filename,'not'), continue; end
    end
    
   
    
    fprintf('\nDoing %s',name);
    
    metrics.name = name;
    
    % load adj matrix
    meta = load([adj_folder,filename]);
    meta = meta.meta;
    
    % load eeg data
    spike = load([eeg_folder,name,not_spike_text,'_eeg.mat']);
    spike = spike.spike;
    
    
    % get sizes for matrices
    n_f = length(meta.spike(1).adj);
    n_times = size(meta.spike(1).index_windows,1);
    n_spikes = length(meta.spike);
    n_ch = size(meta.spike(1).adj(1).adj,2);
    
    % Initialize
    ns = nan(n_spikes,n_times,n_ch,n_f);

    % Initialize spike count
    s_count = 0;
    
    for s = 1:length(meta.spike)

        s_count = s_count + 1;
        
        
        for f = 1:n_f
            adj_all_t= meta.spike(s).adj(f).adj;
            for tt = 1:size(adj_all_t,1)
                
                % Get adj matrix of interest
                adj = squeeze(adj_all_t(tt,:,:));
                
                %% Calculate metrics  
                % node strength
                ns_temp = strengths_und(adj);  
                ns(s,tt,:,f) = ns_temp;
                 
            end
            
        end
           
    end
    
        
  
    metrics.index_windows = meta.spike(1).index_windows;
    metrics.fs = meta.fs;
    metrics.data = ns;
    
    save([out_folder,name,not_spike_text,'_ns.mat'],'metrics');
    
end

end