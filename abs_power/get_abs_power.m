function get_abs_power(time_window,not_a_spike)

%% Description
%{
This function get the absolute signal power in different time windows
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
sig_dev_folder = [results_folder,'power/manual/',time_text];
biggest_dev_folder = [results_folder,'biggest_dev/'];
seq_folder = [results_folder,'seq_data/'];


listing = dir([eeg_folder,'*',not_a_spike_text,'_eeg.mat']);

if exist(sig_dev_folder,'dir') == 0
    mkdir(sig_dev_folder);
end

all_names = {};

% Loop through patients
for i = 1:length(listing)
    
    filename = listing(i).name;
    name_sp = split(filename,'_');
    name = name_sp{1};
    
    % Do either just spikes or just not spikes
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
    values = spike(1).data;
    nchs = size(values,2);
    fs = spike(1).fs;
    n_windows = ntimes;
    
    % Load manual biggest dev file
    if not_a_spike == 0
        manual_big = load([biggest_dev_folder,name,'_rise.mat']);
        manual_big = manual_big.early;
        if length(manual_big.spike) ~= length(spike)
            error('what');
        end
    end
    
    % Load sequence folder
    if not_a_spike == 0
        seq = load([seq_folder,name,'_seq.mat']);
        seq = seq.seq;
        if length(seq) ~= length(spike)
            error('what');
        end
    end
    
    dev_windows = zeros(length(spike),n_windows);
    dev_windows_auto = zeros(length(spike),n_windows);
    dev_windows_rand = zeros(length(spike),n_windows);
    
    dev_windows_first = nan(length(spike),n_windows);
    dev_windows_other = nan(length(spike),n_windows);
    
    % loop through spikes
    for s = 1:length(spike)
        
        %% Initial processing
        % get eeg data
        data = spike(s).data; % ntimes x nch
        
        % pre process
        data = pre_processing(data,do_car,pre_whiten,do_notch,fs);
        
        if not_a_spike == 0
            
            % biggest dev (automatic measurement)
            biggest_dev = spike(s).biggest_dev;

            % manual biggest dev
            biggest_dev_manual = manual_big.spike(s).dev_ch;
            
                       
            % Get sequence info
            if ~isempty(seq(s).seq)
                % first channel in seq
                first_ch = seq(s).first_ch;
                
                % Other channels in sequence (remove first channel)
                other_seq_chs = seq(s).seq(:,1);
                other_seq_chs(other_seq_chs == first_ch) = [];
            end
        end
        
        %% Get power
        % get baseline (diff for each ch)
        baseline = median(data,1); 

        % get deviation from baseline (for all time points) and square to
        % get power
        dev = (abs((data - repmat(baseline,size(data,1),1))).^2); 
        
        % Restrict to biggest dev channel (or average across all channels
        % if not a spike)
        if not_a_spike == 1
            dev_avg_ch = mean(dev,2);
            dev_avg_ch_auto = mean(dev,2);
            
            % pick a random electrode
            nelecs = size(dev,2);
            r_elec = randi(nelecs);
            dev_rand_ch = dev(:,r_elec);
        else
            dev_avg_ch = dev(:,biggest_dev_manual);
            dev_avg_ch_auto = dev(:,biggest_dev);
            
            % For spikes that have a sequence, take the power in the
            % first channel and that in the other channels
            if ~isempty(seq(s).seq)
                dev_first_ch = dev(:,first_ch);
                dev_other_ch = dev(:,other_seq_chs);
            end
        end
        %}
        
        %% Get time windows
        % The peak should be the very center of each file
        peak = round(size(data,1)/2); % size of data should be same as size of values
        
        index_windows = zeros(n_windows,2);
            
        for tt = 1:n_windows
            index_windows(tt,1) = peak + round(time_window(tt)*fs);
            index_windows(tt,2) = peak + round(time_window(tt)*fs) + round(true_window*fs);
        end
        
        % Fix the first and the last to make sure they don't become
        % negative or beyond the total size
        index_windows(1,1) = max(index_windows(1,1),1);
        index_windows(end,2) = min(index_windows(end,2),size(values,1));
        
        
        %% Get the average power in the time window
        % now, get the average power in each time window for that spike
        for t = 1:size(index_windows,1)
            dev_windows(s,t) = mean(dev_avg_ch(max(1,round(index_windows(t,1)))...
                :min(length(dev_avg_ch),round(index_windows(t,2)))));
            
            dev_windows_auto(s,t) = mean(dev_avg_ch_auto(max(1,round(index_windows(t,1)))...
                :min(length(dev_avg_ch_auto),round(index_windows(t,2)))));
            
            if not_a_spike == 1
                dev_windows_rand = mean(dev_rand_ch(max(1,round(index_windows(t,1)))...
                :min(length(dev_rand_ch),round(index_windows(t,2)))));
            end
            
            % first vs other channels in sequence
            if not_a_spike == 0 && ~isempty(seq(s).seq)
                
                dev_windows_first(s,t) = mean(dev_first_ch(max(1,round(index_windows(t,1)))...
                :min(length(dev_first_ch),round(index_windows(t,2)))));
            
            % need to also take the mean across the channels
                dev_windows_other(s,t) = mean(mean(dev_other_ch(max(1,round(index_windows(t,1)))...
                :min(length(dev_other_ch),round(index_windows(t,2))),:)));
            
                
            end
            
        end
        
        if 0
        figure
        set(gcf,'position',[240 264 1201 534]);
        subplot(2,1,1)
        plot(dev_avg_ch);
        hold on
        for t = 1:size(index_windows,1)
            plot([index_windows(t,1) index_windows(t,1)], get(gca,'ylim'),'k--');
            plot([index_windows(t,2) index_windows(t,2)], get(gca,'ylim'),'k--');
        end
        subplot(2,1,2)
        plot(dev_windows(s,:));
        pause
        close
        end
        

    end
    
    % Add to structure
    sig_dev(pt_idx).dev_windows = dev_windows;
    sig_dev(pt_idx).dev_windows_rand = dev_windows_rand;
    sig_dev(pt_idx).dev_windows_auto = dev_windows_auto;
    if not_a_spike == 0 
        sig_dev(pt_idx).dev_windows_first = dev_windows_first;
        sig_dev(pt_idx).dev_windows_other = dev_windows_other;
    else
        sig_dev(pt_idx).dev_windows_first = [];
        sig_dev(pt_idx).dev_windows_other = [];
    end
    sig_dev(pt_idx).fs = fs;
    sig_dev(pt_idx).time_window = time_window;
    sig_dev(pt_idx).index_windows = index_windows;

end

% Save the structure
save([sig_dev_folder,'sig_dev',not_a_spike_text,'.mat'],'sig_dev')


end