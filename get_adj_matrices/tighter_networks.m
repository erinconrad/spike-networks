function tighter_networks(whichPts)

%% Parameters
merge = 1; % merge with existing?
do_car = 1;
pre_whiten = 1;
time_window = 0.5; %0.5 second time window
n_chunks = 11; % How many chunks left and right of spike window


%% Parameters for individualizing spike window
% The time surrounding the spike peak that I will call the spike window if
% I cannot find the signal drop to baseline soon enough
default_spike_window_times = [-1 1];

% Within How many standard deviations away from median must the pre- and
% post-spike period be to be considered baseline
std_allowed = 1;

% How long must the eeg signal be at that baseline to say it's in the pre-
% or post- spike window
time_at_baseline = 0.5; %0.5 seconds

% How long then to move the spike window once I've found this baseline
% period
buffer = 0.05;

freq_bands = [5 15;... %alpha/theta
    15 25;... %beta
    30 40;... % low gamma
    95 105;... % high gamma
    105 256;... %ultra-high
    0 256;... %broadband    
    ]; 
freq_names = {'alpha_theta','beta','low_gamma',...
    'high_gamma','ultra_high','broadband'};


%% Get file locations, load spike times and pt structure
locations = spike_network_files;
main_folder = locations.main_folder;
results_folder = [main_folder,'results/'];
data_folder = [main_folder,'data/'];
script_folder = locations.script_folder;
addpath(genpath(script_folder));
spike_times_file = [data_folder,'spike_times/times.mat'];
pt_file = [data_folder,'spike_structures/pt.mat'];

times = load(spike_times_file); % will result in a structure called "out"
times = times.out;
pt = load(pt_file); % will create a structure called "pt"
pt = pt.pt;

if isempty(whichPts) == 1
    whichPts = 1:length(times);
end

% Loop through patients
for whichPt = whichPts
    
    %% Prep patient
    
    % Skip it if name is empty
    if isempty(times(whichPt).name) == 1, continue; end
    name = times(whichPt).name;
    fprintf('\nDoing %s\n',name);
    
    
    % Get basic data about patient
    pt_folder = [results_folder,name,'/'];
    data = load([pt_folder,'basic_info.mat']); % returns a structure called data
    data = data.data;
    fs = data.fs;
    
     % output folder
    out_folder = [pt_folder,'adj/'];
    if exist(out_folder,'dir') == 0
        mkdir(out_folder);
    end
    
    % convert ch labels to nice form and decide which to ignore
    ch_labels = data.chLabels(:,1);
    ignore = zeros(length(ch_labels),1);
    
    for i = 1:length(ch_labels)
        ch_labels{i} = ieeg_ch_parser(ch_labels{i});
        for j = 1:length(pt(whichPt).ignore.names)
            if strcmp(pt(whichPt).ignore.names(j),ch_labels{i}) == 1
                ignore(i) = 1;
            end
        end
    end
    
    % Loop through spike files
    for f = 1:10
        
        % Skip it if already done
        if merge == 1
            if exist([out_folder,sprintf('finished_%d.mat',f)],'file') ~=0
                fprintf('Did file %d already, skipping...\n',f);
                continue
            end
            
        end
        
        fprintf('Doing file %d of %d...\n',f,10);
        spike = load([pt_folder,sprintf('spikes_%d.mat',f)]);
        spike = spike.spike;
        
        % Initialize output data
        meta_file = [out_folder,sprintf('adj_%d.mat',f)];
        
        % Load it if it exists to see how much we've already done
        if exist(meta_file,'file') ~= 0
            meta = load(meta_file);
            meta = meta.meta;
            
            % Find first unfinished spike
            start_spike = length(meta.spike) + 1;
                
        else
            start_spike = 1;
            meta.name = name;
            meta.which_file = f;
            meta.file_name = sprintf('spikes_%d.mat',f);
        end
        
        
        % Loop through spikes
        for s = start_spike:length(spike)
            fprintf('Doing spike %d of %d...\n',s,length(spike));
            tic
            if isempty(spike(s).time) == 1, continue; end
            
            
            
            % Grab the appropriate channels
            values = spike(s).values(:,~ignore);
            
            % Get spike chs
            is_sp_ch = strcmp(ch_labels(~ignore),spike(s).label);
            is_seq_ch = ismember(ch_labels(~ignore),spike(s).seq_labels);
           
            
            meta.spike(s).time =spike(s).time;
            meta.spike(s).is_sp_ch = is_sp_ch;
            meta.spike(s).is_seq_ch = is_seq_ch;
 
            %% Find nans
            value_nan = isnan(values(:,is_sp_ch));
            values(value_nan,:) = 0;
            
            
            %% Pre-processing
            % Parameters 2 and 3 indicate whether to do CAR and pre-whiten,
            % respectively
            %fprintf('Doing pre-processing...\n');
            values = pre_processing(values,do_car,pre_whiten);      
            values(value_nan,:) = nan;
            nchs = size(values,2);
            
            
            
            %% Decide the spike window
            % The peak should be the very center of each file
            peak = round(size(values,1)/2);
            
            values_seq_ch = values(:,is_sp_ch);
            baseline = nanmedian(values_seq_ch,1);
            dev = abs(values_seq_ch-baseline); % deviation
            max_dev = baseline + nanstd(values_seq_ch,0,1)*std_allowed;
            
            % Find the first point in the spike window: move to the left of
            % the spike, a point at a time, until I've reached a point
            % where there is a minimum duration in which the signal stays
            % around the baseline
            t = 0;
            while 1
                t = t+1; % Increment index moving back
                
                % Specify period
                period = peak - t - time_at_baseline*fs : peak - t;
                
                % Get deviation from baseline in period
                period_dev = dev(period,:);
                
                % See if any indices or channels where it is outside
                % allowed
                outside_dev = period_dev > max_dev;
                
                % If not outside allowed in any, we found a pre-spike
                % period
                if sum(any(outside_dev,2)) == 0
                    spike_window(1) = peak - t - buffer*fs;
                    break
                end
                
                % If t is too far back
                if t > abs(default_spike_window_times(1)*fs)
                    spike_window(1) = peak - abs(default_spike_window_times(1)*fs);
                    break
                end
            end
            
            % Same thing but for last point in spike window
            t = 0;
            while 1
                t = t+1; % Increment index moving back
                
                % Specify period
                period = peak + t : peak + t+time_at_baseline*fs;
                
                % Get deviation from baseline in period
                period_dev = dev(period,:);
                
                % See if any indices or channels where it is outside
                % allowed
                outside_dev = period_dev > max_dev;
                
                % If not outside allowed in any, we found a post-spike
                % period
                if sum(any(outside_dev,2)) == 0
                    spike_window(2) = peak + t + buffer*fs;
                    break
                end
                
                % If t is too far forward
                if t > abs(default_spike_window_times(2)*fs)
                    spike_window(2) = peak + abs(default_spike_window_times(2)*fs);
                    break
                end
            end
            
            
            
            %% Get all the windows
            index_windows = zeros(n_chunks*2+1,2);
            index_windows(n_chunks+1,:) = spike_window; % 12th is spike window
            for i = 1:n_chunks
                
                % when i == 11, this goes from 500 ms before up to the
                % spike window. When i ==1, this goes from 12*500 ms up to
                % 11*500 ms before the spike window
                index_windows(i,:) = [spike_window(1) - time_window*fs*(n_chunks+1-i),...
                    spike_window(1) - time_window*fs*(n_chunks-i)];
                
                % when i == 1, this goes from the spike window to 500 ms
                % after the spike window. When i == 11,...
                index_windows(i+n_chunks+1,:) = [spike_window(2) + time_window*fs*(i-1)...
                    spike_window(2) + time_window*fs*(i)];
            end
            
            % Show the windows
            if 1
                sp_data = values(:,is_sp_ch);
                figure
                plot(linspace(0,14,length(sp_data)),sp_data);
                hold on
                for j = 1:size(index_windows,1)
                    plot([index_windows(j,1) index_windows(j,1)]/fs,get(gca,'ylim'),'k--')
                    plot([index_windows(j,2) index_windows(j,2)]/fs,get(gca,'ylim'),'k--')
                end
                plot([peak peak]/fs,get(gca,'ylim'),'b--');
                plot(get(gca,'xlim'),[baseline baseline],'b--')
                plot(get(gca,'xlim'),[baseline+max_dev baseline+max_dev],'r--');
                plot(get(gca,'xlim'),[baseline-max_dev baseline-max_dev],'r--');
                beep
                pause
                close(gcf)
            end
            
            %% Get adjacency matrices
            %fprintf('Calculating functional networks...\n');

            % Initialize adjacency matrix
            for ff = 1:size(freq_bands,1)
                adj(ff).name = freq_names{ff};
                adj(ff).adj = zeros(n_chunks,nchs,nchs);
            end
            
            for tt = 1:n_chunks*2+1
                
                % get appropriate points
                temp_values = values(round(index_windows(tt,1)):round(index_windows(tt,2)),:); 
                
                % Get adjacency matrices
                t_adj = get_adj_matrices(temp_values,data.fs,freq_bands);
                
                for ff = 1:size(freq_bands,1)
                    adj(ff).adj(tt,:,:) = t_adj(ff).adj;
                end
                
            end
            
            meta.spike(s).adj = adj;
            meta.spike(s).index_windows = index_windows;
            
            % Save the meta file after each spike run
            save(meta_file,'meta');

            t = toc;
            fprintf('Spike %d took %1.1f minutes.\n\n',s,t/60);
            
            if 0
                figure
                set(gcf,'position',[200 250 1175 400]);
                for tt = 1:n_chunks-1
                    subplot(2,5,tt);
                    imagesc(squeeze(adj(4).adj(tt,:,:)));
                    colorbar
                    title(sprintf('%d s',tt))
                end
                beep
                pause
                close(gcf)
            end
        end
        
        % Once I've done all spikes for the file, make a finished.mat
        finished = 1;
        save([out_folder,sprintf('finished_%d.mat',f)],'finished');
        
    end
    
end


end