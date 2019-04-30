function avg_adj_mat(whichPts)

%% Parameters
% 1 = alpha/theta; 2 = beta, 3 = low gamma, 4 = high gamma, 5 = ultra high, 6 = broadband
freq_text = {'alpha/theta','beta','low\ngamma','high\ngamma','ultra high\ngamma','broadband'};
%freq_text = {'alpha/theta'};
n_f = length(freq_text);
n_times = 11;
spike_window_times = [-0.2 0.8];

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
    
    %% Prep avg adjacency matrices
    % Load the first adjacency matrix to get the size in order to
    % initialize one of the arrays
    meta = load([adj_folder,listing(1).name]); 
    meta = meta.meta;
    nchs = length(meta.spike(1).is_seq_ch);
    nfreq = length(meta.spike(1).adj);
    
    for i = 1:nfreq
        adj_avg(i).adj = zeros(n_times,nchs,nchs);
        adj_avg_tiny(i).adj = zeros(n_times,5,5);
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

            s_count = s_count + 1;
            
            seq_chs = meta.spike(s).is_seq_ch; % binary array
            % Just take first five
            non_zero = find(seq_chs);
            if length(non_zero) < 5, continue; end
            non_zero = non_zero(1:5);
            
            seq_chs_5 = ismember(1:length(seq_chs),non_zero);
            
            
            for which_freq = 1:length(meta.spike(s).adj)
                adj_all_t= meta.spike(s).adj(which_freq).adj;  
                adj_avg(which_freq).adj = adj_avg(which_freq).adj + adj_all_t;
                adj_avg_tiny(which_freq).adj = adj_avg_tiny(which_freq).adj + ...
                    adj_all_t(:,seq_chs_5,seq_chs_5);
            end
            
        end
        
    end
    
    %% Divide by number of spikes to get average and calculate global metrics
    ge = zeros(nfreq,n_times);
    sync = zeros(nfreq,n_times);
    
    ge_tiny = zeros(nfreq,n_times);
    sync_tiny = zeros(nfreq,n_times);
    
    for which_freq = 1:length(adj_avg)
        adj_avg(which_freq).adj = adj_avg(which_freq).adj/s_count;
        adj_avg_tiny(which_freq).adj = adj_avg_tiny(which_freq).adj/s_count;
        
        for tt = 1:size(adj_avg(which_freq).adj,1)
            ge(which_freq,tt) = efficiency_wei(squeeze(adj_avg(which_freq).adj(tt,:,:)),0);
            sync(which_freq,tt) = synchronizability_sp(squeeze(adj_avg(which_freq).adj(tt,:,:)));
            
            ge_tiny(which_freq,tt) = efficiency_wei(squeeze(adj_avg_tiny(which_freq).adj(tt,:,:)),0);
            sync_tiny(which_freq,tt) = synchronizability_sp(squeeze(adj_avg_tiny(which_freq).adj(tt,:,:)));
        end
    end
    
    plot_thing(1,:,:) = ge;
    plot_thing(2,:,:) = sync;
    
    plot_title = {'Global efficiency','Synchronizability'};
    
    
    %% Plot aggregated metrics
    
    
    figure
    set(gcf,'position',[1 200 1440 530]);
    [ha, pos] = tight_subplot(n_f-2, n_times, [0 0], [0.08 0.08], [0.05 0.01]);
    for f = 1:n_f-2
        
        % get cmap range within that frequency band
        c_max = max(max(max(adj_avg(f).adj)));
        c_min = min(min(min(adj_avg(f).adj)));
        
        for i = 1:n_times
            axes(ha((f-1)*n_times+i))
            imagesc(squeeze(adj_avg(f).adj(i,:,:)))
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