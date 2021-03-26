function out = get_spike_timing

%% Get file locations, load spike times and pt structure
locations = spike_network_files;
main_folder = locations.main_folder;
results_folder = [main_folder,'results/'];
eeg_folder = [results_folder,'eeg_data/'];

pt_listing = dir([eeg_folder,'*.mat']);
%% Get patient names
all_names = {};
for i = 1:length(pt_listing)
    fname = pt_listing(i).name;
    pt_name = strsplit(fname,'_');
    pt_name = pt_name{1};
    if ~ismember(pt_name,all_names)
        all_names = [all_names;pt_name];
    end
end

% Loop through patients
for p = 1:length(all_names)
    name = all_names{p};
    
    % Load file
    info = load([eeg_folder,name,'_eeg.mat']);
    
    n_spikes = length(info.spike);
    times = zeros(n_spikes,1);
    
    for s = 1:n_spikes
        times(s) = mean([info.spike(s).times(1) info.spike(s).times(2)]);
    end
    
    out(p).name = name;
    out(p).times = times;
end



end