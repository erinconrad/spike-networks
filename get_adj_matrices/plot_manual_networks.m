function plot_manual_networks(simple,brief)

%% Get file locations, load spike times and pt structure
locations = spike_network_files;
main_folder = locations.main_folder;
eeg_folder = [main_folder,'results/eeg_data/'];
results_folder = [main_folder,'results/adj_mat/manual/'];
data_folder = [main_folder,'data/'];
script_folder = locations.script_folder;
addpath(genpath(script_folder));
pt_file = [data_folder,'spike_structures/pt.mat'];


if simple == 1
    out_folder = [results_folder,'adj_simple/'];
    freq_text = {'pairwise\ncorrelation'};
else
    out_folder = [results_folder,'adj_coherence/'];
    freq_text = {'alpha/theta','beta','low\ngamma','high\ngamma'};
end

listing = dir([out_folder,'*_adj.mat']);
new_out_folder = [out_folder,'plots/'];
if exist(new_out_folder,'dir') == 0
    mkdir(new_out_folder);
end

% Loop through finished patients
for i = 1:length(listing)
    name = listing(i).name;
    name_sp = split(name,'_');
    ptname = name_sp{1};
    
    % load
    meta = load([out_folder,name]);
    meta = meta.meta;
    
    % Get number of frequencies
    nfreq = length(meta.spike(1).adj);
    
    % Initialize average networks
    for f = 1:nfreq
        adj_avg(f).adj = zeros(size(meta.spike(1).adj(1).adj));
    end
    
    nspikes = 0;
    % Loop through spikes and contribute to average
    for s = 1:length(meta.spike)
        
        % Skip if sums to zero or if any nans
        if sum(sum(sum(isnan(meta.spike(s).adj(1).adj)))) > 0 || ...
                sum(sum(sum(meta.spike(s).adj(1).adj))) == 0
            continue
        end
        
        nspikes = nspikes + 1;
        
        for f = 1:nfreq
            
            % add current adjaceny to total
            adj_temp = meta.spike(s).adj(f).adj;
            adj_avg(f).adj = adj_avg(f).adj + adj_temp;
               
        end
    end
    
    % Divide by nspikes to get average
    for f = 1:nfreq 
        adj_avg(f).adj = adj_avg(f).adj/nspikes;
    end
    
    % Decide times to loop through
    ntimes = size(adj_avg(1).adj,1);
    if brief == 1
        times = [floor(ntimes/2)-3:floor(ntimes/2)+3];
    else
        times = 1:size(adj_avg(1).adj,1);
    end
    
    % Make plots
    figure
    set(gcf,'position',[1 200 1440 (nfreq-1)*180+180]);
    [ha, pos] = tight_subplot(nfreq, length(times), [0 0], [0.03 0.12], [0.05 0.01]);
    for f = 1:nfreq
        for t =1:length(times)
            axes(ha((f-1)*length(times)+t))
            imagesc(squeeze(adj_avg(f).adj(times(t),:,:)))
            
            if t == 1
                ylabel(sprintf(freq_text{f}));
            end
            
            if f==1
                title(sprintf('%d',t));  
            end
            
            set(gca,'fontsize',15)
            xticklabels([])
            yticklabels([])
        end
    end
    
    print([new_out_folder,ptname],'-depsc');
    close(gcf)
    
end

end