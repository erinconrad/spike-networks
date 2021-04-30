

%% Get file locations, load spike times and pt structure
locations = spike_network_files;
main_folder = locations.main_folder;
hank_data_folder = locations.hank;

%% Load hank data
ev = load([hank_data_folder,'temporalGeneral.mat']);
dat = load([hank_data_folder,'allExpThElecData.mat']);
secs = load([hand_data_folder,'allExpTimingSec.mat']);

npts = length(ev.evStart);

% Loop over cats
for p = 1:npts
    
    %% Get spikes times
    
    % find indices of events with only one spike
    idx = find(ev.numSpikes{p} == 1)';
    
    start = ev.spikeStart{p}'; 
    stop = ev.spikeStop{p}';
    
    % Get start time of these spikes; take this as the peak time (it's
    % slightly earlier)
    spikes = cell2mat(start(idx));
    
    % Data
    values = dat.allExpThElecData{p}';
    
    % Loop through spikes
    for s = 1:length(spikes)
       
        idx = spikes(s);
        
        % Get data clip from 2s before to 0.4 seconds after
        clip_times = [-2 0.4];
        clip_indices = round(clip_times*fs)+idx;
        clip = values(clip_indices(1):clip_indices(2));
        
        % plot the clip
        if 1
            figure
            plot(clip)
            pause
            hold off
        end
        
    end
      
end
