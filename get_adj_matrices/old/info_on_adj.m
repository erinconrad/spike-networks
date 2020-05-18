function info_on_adj(whichPts)

%{
This function just gathers info on what percentage of spikes default to the
2 s spike window and what % use a custom window based on return to baseline
%}


%% Get file locations, load spike times and pt structure
locations = spike_network_files;
main_folder = locations.main_folder;
results_folder = [main_folder,'results/'];
data_folder = [main_folder,'data/'];
script_folder = locations.script_folder;
addpath(genpath(script_folder));
spike_times_file = [data_folder,'spike_times/times.mat'];

times = load(spike_times_file); % will result in a structure called "out"
times = times.out;

if isempty(whichPts) == 1
    whichPts = 1:length(times);
end

for whichPt = whichPts
    
    %% Get basic data about patient
    % Skip it if name is empty
    if isempty(times(whichPt).name) == 1, continue; end
    name = times(whichPt).name;
    fprintf('\nDoing %s\n',name);
    
    pt_folder = [results_folder,name,'/'];
    basic_file = [pt_folder,'basic_info.mat'];
    
    if exist(basic_file,'file') == 0
        continue
    end
    
    basic = load(basic_file);
    fs = basic.data.fs;
    
    adj_folder = [pt_folder,'adj/'];
    
    if exist(adj_folder,'dir') == 0
        continue
    end
    
    sp_windows = nan(1000,2);
    count = 0;
    
    %% Loop through adjacency files
    for f = 1:10
        
        if exist([adj_folder,'adj_',sprintf('%d.mat',f)],'file') == 0
            continue
        end
        
        % Load the file
        adj_file = [adj_folder,'adj_',sprintf('%d.mat',f)];
        adj = load(adj_file);
        meta = adj.meta;
        
        % Loop through the spikes
        for s = 1:length(meta.spike)
            
            count = count + 1;
            % Get the spike window
            sp_window = meta.spike(s).index_windows(12,:)/fs; % 12 is the spike window
            
            sp_windows(count,:) = sp_window;
            
        end
        
    end
    
    
    %% Report stats on spike window for that patient
    % Number of spikes with default spike window
    num_default = sum(abs(diff(sp_windows,1,2) - 2) <= 2e-3); % default spike window is 2 s
    
    mean_window = nanmean(diff(sp_windows,1,2));
    
    median_window = nanmedian(diff(sp_windows,1,2));
    
    fprintf(['For %s, %d spikes out of %d used the default spike window.\n'...
        'The mean and median window duration were %1.3fs and %1.3fs, respectively.\n\n'],...
        name,num_default,count,mean_window,median_window);
    
end

end