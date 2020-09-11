function manual_sig_deviation(time_window,not_a_spike)

%{
This function determines the time periods in which the EEG data surrounding
the spike is significantly different from the baseline (first half second).
This serves as a control for other measures
%}

%% Parameters
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

if length(time_window) == 1
    time_text = sprintf('%1.1f/',time_window);
else
    ntimes = length(time_window);
    all_times = time_window;
    true_window = all_times(2)-all_times(1);
    time_text = sprintf('%1.1f/',true_window);
end

if not_a_spike
    not_a_spike_text = '_not_spike';
else
    not_a_spike_text = '';
end

% Folders
eeg_folder = [results_folder,'eeg_data/'];
sig_dev_folder = [results_folder,'signal_deviation/manual/',time_text];
adj_folder = [results_folder,'adj_mat/manual/adj_simple/',time_text];

listing = dir([eeg_folder,'*',not_a_spike_text,'_eeg.mat']);

if exist(sig_dev_folder,'dir') == 0
    mkdir(sig_dev_folder);
end

all_names = {};
     

for i = 1:length(listing)
    
    filename = listing(i).name;
    name_sp = split(filename,'_');
    name = name_sp{1};
    
    if not_a_spike == 1
        if contains(filename,'not') == 0, continue; end
    else
        if contains(filename,'not') == 1, continue; end
    end
    
    
    [a,b] = ismember(name,all_names);
    if a == 1
        pt_idx = b;
    else
        all_names = [all_names;name];
        pt_idx = length(all_names);
    end

    sig_dev(pt_idx).name = name;
    
    
    fprintf('\nDoing %s...\n',name);
    
    % load eeg data
    spike = load([eeg_folder,filename]);
    spike = spike.spike;
    surround_time = spike(1).surround_time;
    n_spikes = length(spike);
    values = spike(1).data;
    nchs = size(values,2);
    fs = spike(1).fs;
    
    if length(time_window) == 1
        n_windows = round(size(values,1)/fs/time_window); % old way, same total time
    else
        n_windows = ntimes;
    end
    
    dev_windows = zeros(length(spike),n_windows);
    
    % loop through spikes
    for s = 1:length(spike)
        
        % get eeg data
        data = spike(s).data; % ntimes x nch
        
        % pre process
        data = pre_processing(data,do_car,pre_whiten,do_notch,fs);
        
        % get involved chs
        is_sp_ch = spike(s).involved;
        
        % biggest dev
        biggest_dev = spike(s).biggest_dev;
        
        % restrict to involved channels
        %data_spike = data(:,is_sp_ch);
        data_spike = data(:,biggest_dev);
        
        % get baseline (diff for each ch)
        baseline = median(data_spike,1); %1 x n_sp_ch (median across all time points)
        
        % get deviation from baseline (for all time points) and square to
        % get power
        dev = (abs((data_spike - repmat(baseline,size(data_spike,1),1))).^2); % ntimes x n_sp_ch 

        
        % get the average deviation across involved channels
        %dev_avg_ch = mean(dev,2); % ntimes x 1
        dev_avg_ch = dev;
        
        % The peak should be the very center of each file
        peak = round(size(values,1)/2);

        % Get index windows
        if length(time_window) == 1
            index_windows = zeros(n_windows,2);
            tick_window = time_window*fs;

            for tt = 1:n_windows
                index_windows(tt,1) = peak - tick_window*n_windows/2 + tick_window*(tt-1);
                index_windows(tt,2) = peak - tick_window*n_windows/2 + tick_window*(tt);
            end
        else
            index_windows = zeros(n_windows,2);
            
            for tt = 1:n_windows
                index_windows(tt,1) = peak + round(time_window(tt)*fs);
                index_windows(tt,2) = peak + round(time_window(tt)*fs) + round(true_window*fs);
            end
        end
        
        % Fix the first and the last to make sure they don't become
        % negative or beyond the total size
        index_windows(1,1) = max(index_windows(1,1),1);
        index_windows(end,2) = min(index_windows(end,2),size(values,1));
        
        % now, get the average deviation in each time window for that spike
        for t = 1:size(index_windows,1)
            dev_windows(s,t) = mean(dev_avg_ch(max(1,round(index_windows(t,1)))...
                :min(length(dev_avg_ch),round(index_windows(t,2)))));
        end
        
    end
    
    % Now I can compare the signal deviation across time windows
    for t = 2:size(index_windows,1)
        
        % Do a paired t-test
        [~,p,ci,stats] = ttest(dev_windows(:,1),dev_windows(:,t));
        sig_dev(pt_idx).p(t) = p;
        sig_dev(pt_idx).ci(t,:) = ci;
        sig_dev(pt_idx).stats(t) = stats;
        
        
    end
    sig_dev(pt_idx).dev_windows = dev_windows;
    
    % Also get a normalized z-score to combine across patients
    dev_avg_all_spikes = nanmean(dev_windows,1); % avg across spikes
    z_score_dev = (dev_avg_all_spikes - mean(dev_avg_all_spikes))/std(dev_avg_all_spikes);
    sig_dev(pt_idx).z_score_dev = z_score_dev;
    sig_dev(pt_idx).avg_dev = dev_avg_all_spikes;
    sig_dev(pt_idx).time_window = time_window;
    sig_dev(pt_idx).index_windows = index_windows;
    sig_dev(pt_idx).fs = fs;
    
end

% Save the structure
save([sig_dev_folder,'sig_dev',not_a_spike_text,'.mat'],'sig_dev')

end