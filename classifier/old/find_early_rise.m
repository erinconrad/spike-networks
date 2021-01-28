function find_early_rise

%{
Plan:
1) Define the baseline (median) voltage
2) Get the spike peak voltage (highest channel)
3) Define threshold amplitude to be a % of the difference for a certain
amount of time
3) Loop through involved channels
4) Take moving windows and find time in any involved channel above the
threshold amplitude for the minimum time
5) Call the beginning of that the earliest spike rise
%}

%% Parameters
thresh_amp = 0.2;
baseline_time = 2; % first 2 seconds
min_time = 0.02; % 20 ms
overlap = 0.01;
do_notch = 1; % notch filter?
do_car = 1; % common average reference?
pre_whiten = 0; % remove the AR(1) component for a pre-whitening step?

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

eeg_folder = [results_folder,'eeg_data/'];
spike_rise_folder = [results_folder,'spike_rise/'];


if exist(spike_rise_folder,'dir') == 0
    mkdir(spike_rise_folder);
end

listing = dir([eeg_folder,'*','_eeg.mat']);

all_names = {};
     
% Loop through eeg files (contains both spike and not a spike data)
for i = 1:length(listing)
    
    filename = listing(i).name;
    name_sp = split(filename,'_');
    name = name_sp{1};
    
    % skip not a spike files
    if contains(filename,'not') == 1, continue; end
    
    [a,b] = ismember(name,all_names);
    if a == 1
        pt_idx = b;
    else
        all_names = [all_names;name];
        pt_idx = length(all_names);
    end
    
    rise(pt_idx).name = name;
    
    fprintf('\nDoing %s...\n',name);
    
    % load eeg data
    spike = load([eeg_folder,filename]);
    spike = spike.spike;
    n_spikes = length(spike);
    values = spike(1).data;
    nchs = size(values,2);
    fs = spike(1).fs;
    
    % loop through spikes
    for s = 1:length(spike)
        
        % get eeg data
        data = spike(s).data; % ntimes x nch
        
        % pre process
        data = pre_processing(data,do_car,pre_whiten,do_notch,fs);
        
        % get involved chs
        is_sp_ch = spike(s).involved;
        
        % restrict to involved channels
        data_spike = data(:,is_sp_ch);
        
        % get peak time (middle of the file)
        peak = round(size(data,1)/2);
        
        % initialize early time for the spike to be the peak
        early_time = peak;
        
        % Loop through involved channels
        for ich = 1:size(data_spike,2)
            curr_data = data_spike(:,ich);
            
            % Get peak amplitude of that channel
            peak_amp = curr_data(peak);
            
            % Get baseline amplitude of the channel (median over first 2
            % seconds)
            bl_indices = 1:round(baseline_time*fs);
            bl = median(curr_data(bl_indices));
            
            % Get threshold rise
            thresh_rise = abs(peak_amp-bl)*thresh_amp;
            
            % Take moving windows equal to the min time needed to be at the
            % threshold amplitude
            overlap_idx = round(overlap*fs);
            min_idx = round(min_time*fs);
            n_windows = peak/overlap_idx;
            
            v = (0:n_windows)';
            windows = [1+v*overlap_idx,1+v*overlap_idx+min_idx];
            
            early_rise_index = peak; % if find no other time, make it the peak
            % Loop over windows
            for w = 1:size(windows,1)
                amps = curr_data(windows(w,1):windows(w,2));
                diff_from_bl = abs(amps-bl);
                indices_meeting_thresh = diff_from_bl > thresh_rise;
                if sum(indices_meeting_thresh) == length(indices_meeting_thresh)
                    early_rise_index = windows(w,1);
                    break
                end
            end
            
            if early_rise_index < early_time
                early_time = early_rise_index;
            end
            
        end
        
        % plot the spikes and the early time
        if 1
        figure
        offset = 0;
        for ich = 1:size(data_spike,2)
            plot(data_spike(:,ich)+offset);
            if ich<size(data_spike,2)
                offset = offset + max(data_spike(:,ich)) - min(data_spike(:,ich+1));
            end
            hold on
        end
        yl = get(gca,'ylim');
        plot([early_time early_time], yl,'k');
        pause
        close(gcf)
        end
        
        
    end
    
    
end

end