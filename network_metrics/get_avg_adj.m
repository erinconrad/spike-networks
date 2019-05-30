function get_avg_adj(whichPts,small)

%% Parameters
% 1 = alpha/theta; 2 = beta, 3 = low gamma, 4 = high gamma, 5 = ultra high, 6 = broadband
freq_text = {'alpha/theta','beta','low\ngamma','high\ngamma','ultra high\ngamma','broadband'};
%freq_text = {'alpha/theta'};
n_f = length(freq_text);

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
plot_folder = [results_folder,'plots/network_dev/'];
if exist(plot_folder,'dir') == 0
    mkdir(plot_folder);
end

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
    if small == 1
        adj_folder = [results_folder,name,'/adj_small/'];
    else
        adj_folder = [results_folder,name,'/adj/'];
    end
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
    
    n_times = size(index_windows,1);
    
    % Prep adj
    adj_avg = zeros(freq_text,n_times,nchs,nchs);
    
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
            
            for which_freq = 1:length(meta.spike(s).adj)
                adj_all_t= meta.spike(s).adj(which_freq).adj;
                
                for tt = 1:size(adj_all_t,1)
                    
                    % Get adj matrix of interest
                    adj = squeeze(adj_all_t(tt,:,:));
                    
                    if sum(sum(isnan(adj))) > 0
                        continue
                    end
                    
                    adj_avg(which_freq,tt,:,:) = nansum(...
                        squeeze(adj_avg(which_freq,tt,:,:))
                    
                end
                
            end
            
        end
    
    
    
end


end