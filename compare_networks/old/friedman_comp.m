function friedman_comp(whichPts,small)

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
    adj_folder = [results_folder,name,'/adj/'];
    fs = pt(whichPt).fs;
    
    if small == 1
        adj_folder = [results_folder,name,'/adj_small/'];
        stats_folder = [pt_folder,'stats_small/'];
    elseif small == 0
        adj_folder = [results_folder,name,'/adj/'];
        stats_folder = [pt_folder,'stats/'];
    elseif small == 2
        adj_folder = [results_folder,name,'/adj_test/'];
        stats_folder = [pt_folder,'stats_test/'];
    elseif small == 3
        adj_folder = [results_folder,name,'/adj_simple/'];
        stats_folder = [pt_folder,'stats_simple/'];
    elseif small == 4
        adj_folder = [results_folder,name,'/adj_coherence/'];
        stats_folder = [pt_folder,'stats_coherence/'];
    end
    
    
    if exist(stats_folder,'dir') == 0
        mkdir(stats_folder);
    end
    
    %% Load adjacency matrices and calculate metrics
    listing = dir([adj_folder,'adj*.mat']);
    if length(listing) == 0
        fprintf('No adjacency matrices for %s, skipping...\n',name);
        continue
    end
    
    %% Prep avg adjacency matrices
    % Load the first adjacency matrix to get the size in order to
    % initialize one of the arrays
    meta = load([adj_folder,listing(1).name]); 
    meta = meta.meta;
    nchs = length(meta.spike(1).is_seq_ch);
    if small == 3
        nfreq = 1;
        n_times = size(meta.spike(1).adj,1);
    else
        nfreq = length(meta.spike(1).adj);
        n_times = size(meta.spike(1).adj(1).adj,1);
    end
    
    
    for i = 1:nfreq
        adj_avg(i).adj = zeros(n_times,nchs*(nchs-1)/2,1000);
        sim(i).p = zeros(n_times,1);
        sim(i).F = zeros(n_times,1);
  
        
        % First second is the control
        sim(i).p(1) = nan;
        sim(i).F(1) = nan;
    end
    
    
    % Initialize spike count
    s_count = 0;
    
    %% Fill 
    for f = 1:length(listing)
        fname = listing(f).name;
        meta = load([adj_folder,fname]); %outputs a struct named meta
        
        % Load the adjacency matrix file
        meta = meta.meta;
        
        
        % Loop through spikes
        for s = 1:length(meta.spike)

            s_count = s_count + 1;
            for which_freq = 1:nfreq
                if small == 3
                    adj_all_t= meta.spike(s).adj; 
                else
                    adj_all_t= meta.spike(s).adj(which_freq).adj; 
                end
                
                for t = 1:size(adj_all_t,1)
                    % flatten matrix
                    adj_out = flatten_or_expand_adj(squeeze(adj_all_t(t,:,:)));

                    % add to matrix
                    adj_avg(which_freq).adj(t,:,s_count) = adj_out;
                    
                end
               
            end
  
        end
        
    end
    
    %% Remove columns with nans or all zeros
    % Need to add
    for which_freq = 1:nfreq
        adj_avg(which_freq).adj(:,:,any(isnan(adj_avg(which_freq).adj),[1 2])) = [];
    end
    
    for which_freq = 1:nfreq
        adj_avg(which_freq).adj(:,:,sum(adj_avg(which_freq).adj,3)==0) = [];
    end
    
    
    for which_freq = 1:nfreq
        
        % Get the first second (this will be an nch*(nch-1)/2 x 1000
        % matrix)
        first_sec = squeeze(adj_avg(which_freq).adj(1,:,:));
        
        % reshape it for manova (so now it's 1000 x nch*(nch-1)/2)
        first_sec = first_sec';
        
        % make a 1000 x 1 array with time 1
        first_sec_time = ones(size(first_sec,1),1);
        
        % Loop through subsequent times
        for i_times = 2:n_times
            
            fprintf('Doing time %d of %d, freq %d of %d\n',i_times,n_times,which_freq,nfreq);
            
            % Get current time
            test_sec = squeeze(adj_avg(which_freq).adj(i_times,:,:));
                
            % Reshape it for MANOVA
            test_sec = test_sec';
            next_sec_time = 2*ones(size(test_sec,1),1);
            
            % Do perMANOVA
            dis = f_dis([first_sec;test_sec]);
            result = f_permanova(dis,[first_sec_time;next_sec_time],1e3);
            
            sim(which_freq).p(i_times) = result.p;
            sim(which_freq).F(i_times) = result.F;
            
        end
        
        
    end
    
     %% Save the sim structure
    save([stats_folder,'sim_permanova.mat'],'sim');
    
end

end