function plot_avg_adj(whichPts,small)

%% Parameters
times_to_plot = [9:14];

% 1 = alpha/theta; 2 = beta, 3 = low gamma, 4 = high gamma, 5 = ultra high, 6 = broadband
freq_text = {'alpha/theta','beta','low\ngamma','high\ngamma'};
%freq_text = {'alpha/theta'};
if small == 3
    n_f = 1;
else
    n_f = length(freq_text);
end
if small == 3 || small == 4
    n_times = 11*2;
else
    n_times = 11*2 + 1;
end


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
    plot_folder = [results_folder,name,'/plot/'];
    if exist('plot_folder','dir') == 0
        mkdir(plot_folder)
    end
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
        nfreq = n_f;
    end
    
    for i = 1:nfreq
        adj_avg(i).adj = zeros(n_times,nchs,nchs);
    end
    
    % Initialize spike count
    s_count = 0;
    
    for f = 1:length(listing)
        fname = listing(f).name;
        meta = load([adj_folder,fname]); %outputs a struct named meta
        
        % Load the adjacency matrix file
        meta = meta.meta;
        
        
        % Loop through spikes
        for s = 1:length(meta.spike)


            for which_freq = 1:nfreq
                if small == 3
                    if sum(sum(sum(isnan(meta.spike(s).adj)))) > 0
                        continue
                    end
                    adj_all_t= meta.spike(s).adj; 
                else
                    if sum(sum(sum(isnan(meta.spike(s).adj(which_freq).adj)))) > 0
                        continue
                    end
                    adj_all_t= meta.spike(s).adj(which_freq).adj; 
                end
                
                adj_avg(which_freq).adj = adj_avg(which_freq).adj + adj_all_t;
            end
            
            s_count = s_count + 1;
            
        end
        
    end
    
    %% Divide by number of spikes to get average
   
    for which_freq = 1:length(adj_avg)
        adj_avg(which_freq).adj = adj_avg(which_freq).adj/s_count;
        
    end
    
    
    %% Plot average adjacency matrix
    fprintf('Spike count is %d\n',s_count);
    
    
    figure
    set(gcf,'position',[1 200 1440 (nfreq-1)*250+300]);
    [ha, pos] = tight_subplot(nfreq, length(times_to_plot), [0 0], [0.03 0.1], [0.05 0.01]);
    for f = 1:nfreq
        
        % get cmap range within that frequency band
        c_max = max(max(max(adj_avg(f).adj-adj_avg(f).adj(1,:,:))));
        c_min = min(min(min(adj_avg(f).adj-adj_avg(f).adj(1,:,:))));
        
        for i = 1:length(times_to_plot)
            axes(ha((f-1)*length(times_to_plot)+i))
            imagesc(squeeze(adj_avg(f).adj(times_to_plot(i),:,:)))
            %caxis([c_min c_max])
            if i == 1 && small ~= 3
                ylabel(sprintf(freq_text{f}));
            end
            
            if f==1
                title(sprintf('Time window %d',times_to_plot(i)));  
            end
            
            set(gca,'fontsize',20)
            xticklabels([])
            yticklabels([])
        end
    end
    if small == 3
        filename = [name,'_avg_adj_simple'];
    elseif small == 4
        filename = [name,'_avg_adj_coherence'];
    end
    
    print([plot_folder,filename],'-depsc');
    close(gcf)
  
    

end

end