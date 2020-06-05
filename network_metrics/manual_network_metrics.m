function manual_network_metrics(overwrite,simple,time_window)

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
time_text = sprintf('%1.1f/',time_window);


% EEG data folder
eeg_folder = [results_folder,'eeg_data/'];

% Adj mat folder
if simple == 1
    adj_folder = [results_folder,'adj_mat/manual/adj_simple/',time_text];
    network_folder = [results_folder,'networks/manual/simple/',time_text];
end

if exist(network_folder,'dir') == 0
    mkdir(network_folder);
end

listing = dir([adj_folder,'*_adj.mat']);

for i = 1:length(listing)
    
    filename = listing(i).name;
    name_sp = split(filename,'_');
    name = name_sp{1};
    
    if overwrite == 0
        if exist([network_folder,name,'_network_stats.mat'],'file') ~= 0
            fprintf('Already did %s, skipping...\n',name);
            continue;
        end
    end
    
    metrics.name = name;
    
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
                ge(which_freq,s_count,tt) = efficiency_wei(adj,0); % global efficiency of full matrix
                
                % node strength of all chs involved
                ns_temp = strengths_und(adj); 
                ns_seq(which_freq,s_count,tt) = mean(ns_temp(involved));
                
            end
            
        end
            
    end
    
    %% Turn them into z-scores
    z_ns = (((ns_seq-mean(ns_seq,3))./std(ns_seq,0,3)));
    z_ge = (((ge-mean(ge,3))./std(ge,0,3)));
    avg_z_ns = squeeze(nanmean(z_ns,2));
    avg_z_ge = squeeze(nanmean(z_ge,2));
    
    metrics.metric(1).name = 'ns of involved channels';
    metrics.metric(2).name = 'ge';
    metrics.metric(1).val = ns_seq;
    metrics.metric(2).val = ge;
    metrics.metric(1).z = z_ns;
    metrics.metric(2).z = z_ge;
    
    
    %% Significance testing
    
    % Loop through the metrics we're testings
    for m = 1:length(metrics.metric)
        
        % alpha is .05 divided by the number of comparisons
        metrics.metric(m).alpha = 0.05/((size(z_ns,3)-1)*size(z_ns,1));
            
        % Loop through times 2:end and compare each to first time
        for t = 2:size(z_ns,3)
   
            % Loop through the frequencies
            for f = 1:size(metrics.metric(m).z,1)
            
                % Do a 2-sample independent t-test to compare the z scores in each
                [~,p,ci,stats] = ttest2(metrics.metric(m).z(f,:,1),...
                    metrics.metric(m).z(f,:,t));
                metrics.metric(m).p(f,t) = p;
                
                metrics.metric(m).ci(f,t,:) = ci;
                metrics.metric(m).stats(f,t) = stats;
                
            end
            
        end
        
        
    end
    
    %% Save the output structure
    save([network_folder,name,'_network_stats.mat'],'metrics');
    
    %% Plot z-scores over time
    if size(z_ns,1) == 1
        figure
        subplot(2,1,1)
        plot(avg_z_ns,'linewidth',2)
        hold on
        plot(find((metrics.metric(1).p<metrics.metric(1).alpha)),...
            avg_z_ns(metrics.metric(1).p<metrics.metric(1).alpha),...
            'r*');
        title('node strength','fontsize',20)

        subplot(2,1,2)
        plot(avg_z_ge,'linewidth',2)
        hold on
        plot(find((metrics.metric(2).p<metrics.metric(2).alpha)),...
            avg_z_ge(metrics.metric(2).p<metrics.metric(2).alpha),...
            'r*');
        title('ge','fontsize',20)
        if exist([network_folder,'plots/'],'dir') == 0
            mkdir([network_folder,'plots/'])
        end
        print(gcf,[network_folder,'plots/',name],'-depsc');       
        close(gcf)
    end
    
end

end