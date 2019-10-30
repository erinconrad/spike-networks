function detect_network_change(whichPts,small)

%{
As opposed to looking for changes in network statistics, what if I can look
for a more general change in the overall network?

I have a bunch of adjacency matrices at each time point. I could obtain
some measure of correlation across all the adjacency matrices and see if
the adjacency matrices at a certain time have a higher intra-group
correlation than inter-group correlation
%}

n_boot = 1e3;

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
    else
        nfreq = length(meta.spike(1).adj);
    end
    n_times = size(meta.spike(1).adj,1);
    
    for i = 1:nfreq
        adj_avg(i).adj = zeros(n_times,nchs*(nchs-1)/2,1000);
        sim(i).true_sim = zeros(n_times,1);
        sim(i).boot_sim = zeros(n_times,n_boot);
        sim(i).p = zeros(n_times,1);
        
        % First second is the control
        sim(i).boot_sim(1,:) = nan;
        sim(i).p(1) = nan;
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
    
    %% Get intra-group similarity
    for which_freq = 1:nfreq
        for i_times = 1:n_times
            sim(which_freq).true_sim(i_times) = ...
                intra_group_similarity(squeeze(adj_avg(which_freq).adj(i_times,:,:)));   
        end   
    end
    
    %% Permutation test to see if each single time period group has more within than without similarity
    % Use the first time as the one I will shuffle with each subsequent
    % second. This is the one I am using for comparison.
    for which_freq = 1:nfreq
        first_sec = squeeze(adj_avg(which_freq).adj(1,:,:));
        n_first_sec = size(first_sec,2);
        
        % Loop through subsequent times
        for i_times = 2:n_times
            
            % Bootstrap loops
            for i_boot = 1:n_boot
                
                if mod(i_boot,100) == 0
                    fprintf('Doing iteration %d of %d for time %d, %s...\n\n',i_boot,n_boot,i_times,name);
                end
                
                % Get current time
                test_sec = squeeze(adj_avg(which_freq).adj(i_times,:,:));
                
                % Combine vectors from time 1 and current time
                all_vecs = [first_sec,test_sec];
                n_all = size(all_vecs,2);
                
                % Pick a random n_first_sec of then
                p = randperm(n_all,n_first_sec);
                
                % These will form my group
                group = all_vecs(:,p);
                
                % Get intra-group sim
                group_sim = intra_group_similarity(group);
                
                % Fill up the bootstrap vector
                sim(which_freq).boot_sim(i_times,i_boot) = group_sim;
                
            end
            
            % Get a one-sided p-value (the percentage of bootstrap
            % permutations with a similarity greater than or equal to the
            % true one)
            pval = (sum(sim(which_freq).boot_sim(i_times,:) >= ...
                sim(which_freq).true_sim(i_times))+1)/...
                (n_boot+1);
            sim(which_freq).p(i_times) = pval;
            
            % Plot the boot strap
            if 0
                figure
                plot(sim(which_freq).boot_sim(i_times,:),'o')
                hold on
                xl = get(gca,'xlim');
                plot([xl(1) xl(2)],[sim(which_freq).true_sim(i_times) sim(which_freq).true_sim(i_times)])
                title(sprintf('Time %d p-value = %1.3f',i_times,pval))
                pause
                close(gcf)
                
            end
        end
    end
    
    %% Save the sim structure
    save([stats_folder,'sim.mat'],'sim');
        
    
end


end

