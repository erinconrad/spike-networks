function spike_pca(whichPts)

n_f = 6;

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
plot_folder = [results_folder,'plots/just_spike/'];

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
    adj_folder = [results_folder,name,'/just_spike/'];
    fs = pt(whichPt).fs;
    
    stats_folder = [pt_folder,'stats/'];
    if exist(stats_folder,'dir') == 0
        mkdir(stats_folder);
    end
    
    %% Load adjacency matrices and calculate metrics
    listing = dir([adj_folder,'adj*.mat']);
    if length(listing) == 0
        fprintf('No adjacency matrices for %s, skipping...\n',name);
        continue
    end
    
    % Number of spikes expected
    n_spikes = sum(~isnan(times(whichPt).spike_times));
    
    % Load the first adjacency matrix to get the size in order to
    % initialize one of the arrays
    meta = load([adj_folder,listing(1).name]); 
    meta = meta.meta;
    nchs = length(meta.spike(1).is_seq_ch);
    
    % Initialize matrices
    sp_adj = nan(n_f,n_spikes,nchs*(nchs-1)/2);
    sp_times = nan(n_spikes,1);
    sp_labels = cell(n_spikes,1);
    
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
            
            
            % Get spike time and ch
            sp_times(s_count) = spike(s_count).time;
            sp_labels{s_count} = spike(s_count).label;
            
            for which_freq = 1:length(meta.spike(s).adj)
                adj_all_t= meta.spike(s).adj(which_freq).adj;
                
                % Get the spikey part
                adj = squeeze(adj_all_t(2,:,:));
                
                % Flatten just the upper triangle
                adj_flat = zeros(nchs*(nchs-1)/2,1);
                flat_idx = 0;
                
                for i = 1:nchs
                    for j = 1:i-1
                        flat_idx = flat_idx + 1;
                        adj_flat(flat_idx) = adj(j,i);
                    end
                end
                
                % add to output matrix
                sp_adj(which_freq,s_count,:) = adj_flat;
                
                
            end
        end
        
    end
    
    % Output to file
    sp_adj_info.adj = sp_adj;
    sp_adj_info.time = sp_times;
    sp_adj_info.labels = sp_labels;
    
    save([stats_folder,'sp_adj.mat'],'sp_adj_info');
    
end


end