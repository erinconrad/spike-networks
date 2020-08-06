function plot_avg_spike

%% Parameters
do_notch = 1; % notch filter?
do_car = 1; % common average reference?
pre_whiten = 0; % remove the AR(1) component for a pre-whitening step?

%% Get file locations, load spike times and pt structure
locations = spike_network_files;
main_folder = locations.main_folder;
results_folder = [main_folder,'results/'];
data_folder = [main_folder,'data/'];
eeg_folder = [main_folder,'results/eeg_data/'];
script_folder = locations.script_folder;
addpath(genpath(script_folder));
pt_file = [data_folder,'spike_structures/pt.mat'];
bct_folder = locations.BCT;
addpath(genpath(bct_folder));
out_folder = [results_folder,'plots/'];
sig_dev_folder = [results_folder,'signal_deviation/manual/'];


if exist(out_folder,'dir') == 0
    mkdir(out_folder);
end

% Load spike file for one patient to get the surround time
spike = load([eeg_folder,'HUP074_eeg.mat']);
spike = spike.spike;
surround_time = spike(1).surround_time;

% get full directory listing
listing = dir(eeg_folder);
count = 0;
pt_names = {};
for i = 1:length(listing)
    
    % Get name
    fname = listing(i).name;
    if contains(fname,'_eeg.mat') == 0, continue; end
    pt_name = strsplit(fname,'_');
    pt_name = pt_name{1};
    pt_names = [pt_names;pt_name];
    
    % Load the file
    spike = load([eeg_folder,fname]);
    spike = spike.spike;
    nspikes = length(spike);
    
    dev_all = zeros(size(spike(1).data,1),nspikes);
    surround = spike(1).surround_time;
    fs = spike(1).fs;
    
    for s = 1:length(spike)
        involved = spike(s).involved;
        
        values = spike(s).data;
        
        % pre process
        values = pre_processing(values,do_car,pre_whiten,do_notch,fs);
        
        x = values(:,involved);
        
        % get baseline
        bl = median(x,1);
        
        % get dev
        dev = abs(x-bl);
        
        % avg across spike chs
        avg_dev = mean(dev,2);
        
        dev_all(:,s) = avg_dev;
        
    end
    
    final_dev = mean(dev_all,2);
    plot(linspace(-surround,surround,length(final_dev)),final_dev)
    pause
        
end


end