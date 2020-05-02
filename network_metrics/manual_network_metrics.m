function manual_network_metrics(simple)

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

% EEG data folder
eeg_folder = [results_folder,'eeg_data/'];

% Adj mat folder
if simple == 1
    adj_folder = [results_folder,'adj_mat/manual/adj_simple/'];
    network_folder = [results_folder,'manual/simple/plots/'];
end

if exist(network_folder,'dir') == 0
    mkdir(network_folder);
end

listing = dir([adj_folder,'*_adj.mat']);

for i = 1:length(listing)
    
    filename = listing(i).name;
    name_sp = split(filename,'_');
    name = name_sp{1};
    
    % load adj matrix
    meta = load([adj_folder,filename]);
    meta = meta.meta;
    
    % load eeg data
    spike = load([eeg_folder,name,'_eeg.mat']);
    spike = spike.spike;
    
    % get sizes for matrices
    n_f = length(meta.spike(1).adj);
    n_times = size(meta.spike(1).index_windows,1);
    n_spikes = length(meta.spike);
    
    % Initialize
    ge = nan(n_f,n_spikes,n_times);
    ns_seq = nan(n_f,n_spikes,n_times);
    
    % Initialize spike count
    s_count = 0;
    out = [];
    
    for s = 1:length(meta.spike)

        s_count = s_count + 1;

        fprintf('Doing spike %d of %d\n',s_count,n_spikes);
        involved = spike(s).involved;
        
        for which_freq = 1:n_f
            adj_all_t= meta.spike(s).adj(which_freq).adj;
            for tt = 1:size(adj_all_t,1)
                
                % Get adj matrix of interest
                adj = squeeze(adj_all_t(tt,:,:));
                
                %% Calculate metrics
                ge(which_freq,s_count,tt) = efficiency_wei(adj,0); % global efficiency
                
                % node strength of all ch involved
                ns_temp = strengths_und(adj); 
                ns_seq(which_freq,s_count,tt) = mean(ns_temp(involved));
                
            end
            
        end
            
    end
    
    z_ns = (((ns_seq-mean(ns_seq,3))./std(ns_seq,0,3)));
    z_ge = (((ge-mean(ge,3))./std(ge,0,3)));
    avg_z_ns = squeeze(nanmean(z_ns,2));
    avg_z_ge = squeeze(nanmean(z_ge,2));
    
    figure
    subplot(2,1,1)
    plot(avg_z_ns)
    title('node strength')
    
    subplot(2,1,2)
    plot(avg_z_ge)
    title('ge')
    print(gcf,[network_folder,name],'-depsc');       
    close(gcf)
    
end

end