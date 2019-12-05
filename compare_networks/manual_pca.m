function manual_pca(whichPts,simple)

%% Parameters
sp_time = 11;

%% Get file locations, load spike times and pt structure
locations = spike_network_files;
main_folder = locations.main_folder;
results_folder = [main_folder,'results/'];
data_folder = [main_folder,'data/'];
script_folder = locations.script_folder;
addpath(genpath(script_folder));
pt_file = [data_folder,'spike_structures/pt.mat'];
bct_folder = locations.BCT;
addpath(genpath(bct_folder));
plot_folder = [results_folder,'plots/'];
if exist(plot_folder,'dir') == 0
    mkdir(plot_folder);
end


pt = load(pt_file); % will create a structure called "pt"
pt = pt.pt;

for whichPt = whichPts

    % Skip it if name is empty
    if isempty(pt(whichPt).name) == 1, continue; end
    name = pt(whichPt).name;
    
    if simple == 1
        adj_folder = [results_folder,'adj_mat/adj_simple/'];
    elseif simple == 0
        adj_folder = [results_folder,'adj_mat/adj_coherence/'];
    end

    
    %% Load adjacency matrices 
    meta = load([adj_folder,name,'_adj.mat']);
    meta = meta.meta;
    
    % Get number of spikes, frequencies, times, channels
    nspikes = length(meta.spike);
    nfreq = length(meta.spike(1).adj);
    n_times = size(meta.spike(1).adj(1).adj,1);
    nchs = size(meta.spike(1).adj(1).adj,2);
    
    %% Get arrays of flattened adjacency matrices
    for i = 1:nfreq
        adj_avg(i).adj = zeros(n_times,nchs*(nchs-1)/2,nspikes);
    end
    
    % Loop over spikes
    for s = 1:length(meta.spike)
        
        for which_freq = 1:nfreq
        
            adj_all_t= meta.spike(s).adj(which_freq).adj;

            for t = 1:size(adj_all_t,1)
                % flatten matrix
                adj_out = flatten_or_expand_adj(squeeze(adj_all_t(t,:,:)));

                % add to matrix
                adj_avg(which_freq).adj(t,:,s) = adj_out;

            end
            
        end
        
    end
    
    %% Remove columns with any nans or all zeros
    for which_freq = 1:nfreq
        adj_avg(which_freq).adj(:,:,any(isnan(adj_avg(which_freq).adj),[1 2])) = [];
    end
    
    for which_freq = 1:nfreq
        adj_avg(which_freq).adj(:,:,sum(adj_avg(which_freq).adj,3)==0) = [];
    end
    
    %% Look at just the spike time and do pca
    for f = 1:nfreq
        
        % Transpose so that it is nspikes x nchannels
        sp_adj(f).adj = squeeze(adj_avg(f).adj(sp_time,:,:))';
        
        % do pca
        [coeff,score,latent] = pca(sp_adj(f).adj);
        
        % Get the first 2 subgraphs
        for i = 1:2
            
            % the flattened adj matrix for the principal component
            subgraph = coeff(:,i);
        
            % expand to get full adjacency matrix
            A = flatten_or_expand_adj(subgraph);
            nchs = size(A,1);
            
            % Plot the histogram of scores of the component for the 1000
            % spikes
            figure
            imagesc(A)
            title(sprintf('Adjacency matrix for component %d',i))
            
        end
        
    end
    
    
end

end