function manual_sig_deviation(time_window)

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
time_text = sprintf('%1.1f/',time_window);

% Folders
eeg_folder = [results_folder,'eeg_data/'];
sig_dev_folder = [results_folder,'signal_deviation/manual/',time_text];
adj_folder = [results_folder,'adj_mat/manual/adj_simple/',time_text];

listing = dir([adj_folder,'*_adj.mat']);

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
    surround_time = spike(1).surround_time;
    
    % load adjacency matrix data (which will allow me to get the index
    % windows
    meta = load([adj_folder,name,'_adj.mat']);
    meta = meta.meta;
    index_windows = meta(1).spike(1).index_windows;
    sig_dev(i).index_windows = index_windows;
    sig_dev(i).surround_time = surround_time;
    
    dev_windows = zeros(length(spike),size(index_windows,1));
    fs = spike(1).fs;
    
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
        data_spike = data(:,is_sp_ch);
        %data_spike = data(:,biggest_dev);
        
        % get baseline (diff for each ch)
        baseline = median(data_spike,1); %1 x n_sp_ch (median across all time points)
        
        % get deviation from baseline (for all time points) and square to
        % get power
        dev = (abs((data_spike - repmat(baseline,size(data_spike,1),1))).^2); % ntimes x n_sp_ch 

        
        % get the average deviation across involved channels
        dev_avg_ch = mean(dev,2); % ntimes x 1
        %dev_avg_ch = dev;
        
        % now, get the average deviation in each time window for that spike
        for t = 1:size(index_windows,1)
            dev_windows(s,t) = mean(dev_avg_ch(round(index_windows(t,1)):round(index_windows(t,2))));
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
    
    % Also get a normalized z-score to combine across patients
    dev_avg_all_spikes = nanmean(dev_windows,1); % avg across spikes
    z_score_dev = (dev_avg_all_spikes - mean(dev_avg_all_spikes))/std(dev_avg_all_spikes);
    sig_dev(i).z_score_dev = z_score_dev;
    sig_dev(i).avg_dev = dev_avg_all_spikes;
    
    % Save the structure
    save([sig_dev_folder,'sig_dev.mat'],'sig_dev')
    
end

end