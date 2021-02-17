function get_ns(overwrite,simple,time_window,not_spike)

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
biggest_dev_folder = [results_folder,'biggest_dev/'];
seq_folder = [results_folder,'seq_data/'];



% Adj mat folder
if simple == 1
    adj_folder = [results_folder,'adj_mat/manual/adj_simple/',time_text];
    metrics_folder = [results_folder,'metrics/manual/simple/',time_text];
else
    adj_folder = [results_folder,'adj_mat/manual/adj_coherence/',time_text];
    metrics_folder = [results_folder,'metrics/manual/coherence/',time_text];
end

if exist(metrics_folder,'dir') == 0
    mkdir(metrics_folder);
end

if not_spike == 1
    not_spike_text = '_not_spike';
else
    not_spike_text = '';
end

listing = dir([adj_folder,'*_adj.mat']);

for i = 1:length(listing)
    
    filename = listing(i).name;
    name_sp = split(filename,'_');
    name = name_sp{1};
    if not_spike
        if ~contains(filename,'not'), continue; end
    else
        if contains(filename,'not'), continue; end
    end
    
    if overwrite == 0
        if exist([metrics_folder,name,not_spike_text,'_ns.mat'],'file') ~= 0
            fprintf('Already did %s, skipping...\n',name);
            continue;
        end
    end
    
    fprintf('\nDoing %s',name);
    
    metrics.name = name;
    
    % load adj matrix
    meta = load([adj_folder,filename]);
    meta = meta.meta;
    
    % load eeg data
    spike = load([eeg_folder,name,not_spike_text,'_eeg.mat']);
    spike = spike.spike;
    
    % Load manual biggest dev file
    if not_spike == 0
        manual_big = load([biggest_dev_folder,name,'_rise.mat']);
        manual_big = manual_big.early;
        if length(manual_big.spike) ~= length(spike)
            error('what');
        end
    end
    
    % Load sequence folder
    if not_spike == 0
        seq = load([seq_folder,name,'_seq.mat']);
        seq = seq.seq;
        if length(seq) ~= length(spike)
            error('what');
        end
    end
    
    % get sizes for matrices
    n_f = length(meta.spike(1).adj);
    n_times = size(meta.spike(1).index_windows,1);
    n_spikes = length(meta.spike);
    n_ch = size(meta.spike(1).adj(1).adj,2);
    
    % Initialize
    ns = nan(n_f,n_spikes,n_times);
    ns_auto = nan(n_f,n_spikes,n_times);
    ns_all = nan(n_f,n_spikes,n_times,n_ch);
    ns_first = nan(n_f,n_spikes,n_times);
    ns_other = nan(n_f,n_spikes,n_times);
    ge = nan(n_f,n_spikes,n_times);
    trans = nan(n_f,n_spikes,n_times);
    involved = nan(n_spikes,n_ch);

    % Initialize spike count
    s_count = 0;
    
    for s = 1:length(meta.spike)

        s_count = s_count + 1;
        
        %fprintf('Doing spike %d of %d\n',s_count,n_spikes);
        if not_spike == 0
            biggest_dev = spike(s).biggest_dev;
            biggest_dev_manual = manual_big.spike(s).dev_ch;
            involved(s,:) = spike(s).involved;
            
            if ~isempty(seq(s).seq)
                % first channel in seq
                first_ch = seq(s).first_ch;
                other_seq_chs = seq(s).seq(:,1);
                other_seq_chs(other_seq_chs == first_ch) = [];
            end
        end
        
        for f = 1:n_f
            adj_all_t= meta.spike(s).adj(f).adj;
            for tt = 1:size(adj_all_t,1)
                
                % Get adj matrix of interest
                adj = squeeze(adj_all_t(tt,:,:));
                
                %% Calculate metrics  
                % node strength of biggest dev channel
                ns_temp = strengths_und(adj);  
                if not_spike == 0
                    ns(f,s,tt) = ns_temp(biggest_dev_manual);
                    ns_auto(f,s,tt) = ns_temp(biggest_dev);
                    if ~isempty(seq(s).seq)
                        ns_first(f,s,tt) = ns_temp(first_ch);
                        ns_other(f,s,tt) = mean(ns_temp(other_seq_chs));
                    end
                else
                    ns(f,s,tt) = mean(ns_temp); % average ns if not spike
                    ns_auto(f,s,tt) = mean(ns_temp);
                    
                end
                ns_all(f,s,tt,:) = ns_temp;
                
                % global efficiency of full matrix
                ge(f,s,tt) = efficiency_wei(adj,0); 
                
                % transitivity of the full matrix
                trans(f,s,tt) = transitivity_wu(adj);
                
            end
            
        end
           
    end
    
        
    % Fill structure
    for f = 1:n_f
        if isfield(meta.spike(1).adj(f),'name') == 1
            metrics.freq(f).name = meta.spike(1).adj(f).name;
        else
            metrics.freq(f).name = 'correlation';
        end
        metrics.freq(f).ns.name = 'node strength biggest channel';
        metrics.freq(f).ns.data = squeeze(ns(f,:,:));
        
        metrics.freq(f).ns_auto.name = 'node strength auto biggest channel';
        metrics.freq(f).ns_auto.data = squeeze(ns_auto(f,:,:));
        
        metrics.freq(f).ge.name = 'global efficiency';
        metrics.freq(f).ge.data = squeeze(ge(f,:,:));
        
        metrics.freq(f).trans.name = 'transitivity';
        metrics.freq(f).trans.data = squeeze(trans(f,:,:));
        
        metrics.freq(f).ns_all.name = 'node strength all';
        metrics.freq(f).ns_all.data = squeeze(ns_all(f,:,:,:));

        metrics.freq(f).involved = involved;
        
        if not_spike == 0
            metrics.freq(f).ns_first.data = squeeze(ns_first(f,:,:));
            metrics.freq(f).ns_other.data = squeeze(ns_other(f,:,:));
        end

    end
    if not_spike == 0
        metrics.biggest_dev = biggest_dev;
    else
        metrics.biggest_dev = nan;
    end
    metrics.index_windows = meta.spike(1).index_windows;
    metrics.fs = meta.fs;
    
    save([metrics_folder,name,not_spike_text,'_ns.mat'],'metrics');
    
end

end