function avg_adj_mat(whichPts,small)

%% Parameters
% 1 = alpha/theta; 2 = beta, 3 = low gamma, 4 = high gamma, 5 = ultra high, 6 = broadband
freq_text = {'alpha/theta','beta','low\ngamma','high\ngamma','ultra high\ngamma','broadband'};
%freq_text = {'alpha/theta'};
n_f = length(freq_text);
n_times = 11*2+1;
do_plot = 0;

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
    nfreq = length(meta.spike(1).adj);
    
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

            if sum(sum(sum(isnan(meta.spike(s).adj(1).adj)))) > 0
                continue
            end
            
       
            
            for which_freq = 1:length(meta.spike(s).adj)
                adj_all_t= meta.spike(s).adj(which_freq).adj;  
                adj_avg(which_freq).adj = adj_avg(which_freq).adj + adj_all_t;
            end
            
            s_count = s_count + 1;
            
        end
        
    end
    
    %% Divide by number of spikes to get average and calculate global metrics
    ge = nan(nfreq,n_times);
    sync = nan(nfreq,n_times);
    
   
    
    for which_freq = 1:length(adj_avg)
        adj_avg(which_freq).adj = adj_avg(which_freq).adj/s_count;
        
        if sum(sum(sum(isnan(adj_avg(which_freq).adj)))) > 0
            continue
        end
        
        for tt = 1:size(adj_avg(which_freq).adj,1)
            ge(which_freq,tt) = efficiency_wei(squeeze(adj_avg(which_freq).adj(tt,:,:)),0);
            sync(which_freq,tt) = synchronizability_sp(squeeze(adj_avg(which_freq).adj(tt,:,:)));
            
        end
    end
    
    %% Save the avg adjacency matrix
    if small == 0
        save([stats_folder,'/avg_adj.mat'],'adj_avg');
    elseif small == 1
        save([stats_folder,'/avg_adj_small.mat'],'adj_avg');
    elseif small == 2
        save([stats_folder,'/avg_adj_test.mat'],'adj_avg');
    end
    
    plot_thing(1,:,:) = ge;
    plot_thing(2,:,:) = sync;
    
    plot_title = {'Global efficiency','Synchronizability'};
    
    
    %% Plot aggregated metrics
    fprintf('Spike count is %d\n',s_count);
    
    if do_plot == 1

    
    figure
    set(gcf,'position',[1 200 1440 530]);
    [ha, pos] = tight_subplot(n_f-2, n_times, [0 0], [0.08 0.08], [0.05 0.01]);
    for f = 1:n_f-2
        
        % get cmap range within that frequency band
        c_max = max(max(max(adj_avg(f).adj-adj_avg(f).adj(1,:,:))));
        c_min = min(min(min(adj_avg(f).adj-adj_avg(f).adj(1,:,:))));
        
        for i = 1:n_times
            axes(ha((f-1)*n_times+i))
            imagesc(squeeze(adj_avg(f).adj(i,:,:)-adj_avg(f).adj(1,:,:)))
            caxis([c_min c_max])
            if i == 1
                ylabel(sprintf(freq_text{f}));
            end
            
            if f==1
                title(sprintf('%d s',i));  
            end
            
            set(gca,'fontsize',20)
            xticklabels([])
            yticklabels([])
        end
    end
    filename = [name,'_avg_adj'];
    print([plot_folder,filename],'-depsc');
    
    
    if 1
    figure
    set(gcf,'position',[26 0 1242 900])
    [ha, pos] = tight_subplot(n_f-2, size(plot_thing,1), [0.04 0.04], [0.08 0.08], [0.05 0.01]);
    for f = 1:n_f-2
        for i = 1:size(plot_thing,1)
            axes(ha((f-1)*size(plot_thing,1)+i))
            plot(squeeze(plot_thing(i,f,:)),'ks-','linewidth',2)
            if f == 1
                title(sprintf(plot_title{i}));
            end
            if f == n_f-2
               xlabel('Time (s)') 
            end
            yticklabels([])
            
            if i == 1
                ylabel(sprintf(freq_text{f}));
            end
            set(gca,'fontsize',20)
        end
    end
    filename = [name,'_avg_adj_metrics'];
    print([plot_folder,filename],'-depsc');
    end
    
    end
    
    
end


end