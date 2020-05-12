function manual_sig_deviation

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

% Folders
eeg_folder = [results_folder,'eeg_data/'];
sig_dev_folder = [results_folder,'signal_deviation/manual/'];
adj_folder = [results_folder,'adj_mat/manual/adj_simple/'];

listing = dir([eeg_folder,'*_eeg.mat']);

if exist(sig_dev_folder,'dir') == 0
    mkdir(sig_dev_folder);
end

for i = 1:length(listing)
    
    filename = listing(i).name;
    name_sp = split(filename,'_');
    name = name_sp{1};
    sig_dev(i).name = name;
    
    % load eeg data
    spike = load([eeg_folder,name,'_eeg.mat']);
    spike = spike.spike;
    
    % load adjacency matrix data (which will allow me to get the index
    % windows
    meta = load([adj_folder,name,'_adj.mat']);
    meta = meta.meta;
    index_windows = meta(1).spike(1).index_windows;
    sig_dev(i).index_windows = index_windows;
    
    dev_windows = zeros(length(spike),size(index_windows,1));
    
    % loop through spikes
    for s = 1:length(spike)
        
        % get eeg data
        data = spike(s).data; % ntimes x nch
        
        % get involved chs
        is_sp_ch = spike(s).involved;
        
        % restrict to involved channels
        data_spike = spike(s).data(:,is_sp_ch);
        
        % get baseline (diff for each ch)
        baseline = median(data_spike,1); %1 x n_sp_ch
        
        % get deviation from baseline
        dev = sqrt((data_spike - repmat(baseline,size(data_spike,1),1)).^2); % ntimes x n_sp_ch 
        
        % get the average deviation across involved channels
        dev_avg_ch = mean(dev,2); % ntimes x 1
        
        % now, get the average deviation in each time window for that spike
        for t = 1:size(index_windows,1)
            dev_windows(s,t) = mean(dev_avg_ch(index_windows(t,1):index_windows(t,2)));
        end
        
    end
    
    % Now I can compare the signal deviation across time windows
    for t = 2:size(index_windows,1)
        
        % Do a two-sample t-test
        [~,p,ci,stats] = ttest2(dev_windows(:,1),dev_windows(:,t));
        sig_dev(i).p(t) = p;
        sig_dev(i).ci(t,:) = ci;
        sig_dev(i).stats(t) = stats;
        
    end
    
    % Save the structure
    save([sig_dev_folder,'sig_dev.mat'],'sig_dev')
    
end

end