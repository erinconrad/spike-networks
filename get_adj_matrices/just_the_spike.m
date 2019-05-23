function just_the_spike(whichPts)

%% Parameters
merge = 1; % merge with existing?
do_car = 1;
pre_whiten = 0;
non_spike_time = 3; % second 3 of record, chosen arbitrarily

freq_bands = [5 15;... %alpha/theta
    15 25;... %beta
    30 40;... % low gamma
    95 105;... % high gamma
    105 256;... %ultra-high
    0 256;... %broadband    
    ]; 
freq_names = {'alpha_theta','beta','low_gamma',...
    'high_gamma','ultra_high','broadband'};


%% Filter parameters
flow = 4; % low pass for slow wave
fhigh = 40; % high pass for spikey bit
butter_order = 6; % butterworth order
min_slow = 50; % how many points long does slow wave need to be
post_spike_buffer = 50; % define slow wave to start 50 points after spikey part ends


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
    out_folder = [pt_folder,'just_spike/'];
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
    
    %% Build the filters
    fn = fs/2; %Nyquist
    
    % High pass filter for spikey bit
    [B_high,A_high] = butter(butter_order,fhigh/fn,'high');
    
    % Low pass filter to get the slow wave
    [B_low,A_low] = butter(butter_order,flow/fn,'low');
    
    
    %% Loop through spikes and get adjacency matrices
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
            nchs = size(values,2);
            peak = round(size(values,1)/2);
           
            
            %% Common average reference
            value_nan = isnan(values(:,is_sp_ch));
            values(value_nan,:) = 0;
            values = pre_processing(values,1,0);      
            values(value_nan,:) = nan;
            
            
            
            %% Find spikey bit     
            % Isolate the first spike channel
            values_sp = values(:,is_sp_ch);
            
            % High pass filter to get the spikey bit, take abs
            hp = abs(filtfilt(B_high,A_high,values_sp));
            
            % Now lp filter this to decide when it rises above baseline
            hp_lp = filtfilt(B_low,A_low,hp);
            hp_bl = nanmedian(hp_lp);
            hp_std = nanstd(hp_lp);
            above_baseline = hp_lp > hp_bl + hp_std;
            
            % Get spikey start
            pre_sp_bin = above_baseline(1:peak);
            pre_sp_below = find(pre_sp_bin==0);
            spikey_start = pre_sp_below(end);
            
            % Get spikey end
            post_sp_bin = above_baseline(peak:end);
            post_sp_below = find(post_sp_bin==0);
            spikey_end = peak + post_sp_below(1);
            
            %% Find slow wave
            
            % Define start of slow wave to be spikey end
            slow_start = spikey_end + post_spike_buffer;
            
            % Find end of slow wave
            % Low pass filter to get the slow wave
            lp = (filtfilt(B_low,A_low,values_sp));
            lp_bl = nanmedian(lp); % get baseline
            
            % Get if above or below baseline, and look for times it
            % switches
            diff_sign_dev_lp = diff(sign(lp - lp_bl));
            
            % Only take times after the spikey part
            after_spikey = diff_sign_dev_lp(slow_start+min_slow:end);
            
            % Find the indices with a NEGATVIE sign change
            sign_change = find(after_spikey < 0);
            slow_end = slow_start + min_slow+ sign_change(1);
            
             % The indices of the start and end index of the spikey bit and
            % slow bit
            spikey = [spikey_start spikey_end];
            slow = [slow_start slow_end];
            
            
            %% Pick an arbitrary other time for baseline
            % This will sometimes have a spike in it, but that's ok.
            non_spike(1) = non_spike_time*fs;
            
            % Make the end index be this plus the average spikey and slow
            % duration
            non_spike(2) = non_spike(1) + ((spikey(2)-spikey(1))+slow(2)-slow(1))/2;
            
            windows = [non_spike;spikey;slow];
            
            if spikey(2)-spikey(1) < 100
                fprintf('Skipping as spike less than 100 points long\n');
                continue;
            end
            
            if 0
                figure
                set(gcf,'position',[169 548 1272 250])
                plot(values_sp,'linewidth',2)
                hold on
                plot(lp)
                plot(get(gca,'xlim'),[lp_bl lp_bl],'k--')
                
                %plot(hp)
                %plot(hp_lp)
                %plot(above_baseline*max(values_sp))
                for i = 1:size(windows,1)
                    plot([windows(i,1) windows(i,1)],get(gca,'ylim'),...
                        'linewidth',2)
                    plot([windows(i,2) windows(i,2)],get(gca,'ylim'),...
                        'linewidth',2)
                end
                pause
                close(gcf)
                continue
            end
            
   
            
            %% Pre-whiten
            value_nan = isnan(values(:,is_sp_ch));
            values(value_nan,:) = 0;
            values = pre_processing(values,0,1);      
            values(value_nan,:) = nan;
            
            %% Get adjacency matrices for spikey bit, slow bit, and non spike bit
            
            
            % Initialize adjacency matrix
            for ff = 1:size(freq_bands,1)
                adj(ff).name = freq_names{ff};
                adj(ff).adj = zeros(size(windows,1),nchs,nchs);
            end
            
            for tt = 1:size(windows,1)
                
                % get appropriate points
                temp_values = values(round(windows(tt,1)):round(windows(tt,2)),:); 
                
                % Get adjacency matrices
                t_adj = get_adj_matrices(temp_values,fs,freq_bands);
                
                for ff = 1:size(freq_bands,1)
                    adj(ff).adj(tt,:,:) = t_adj(ff).adj;
                end
                
            end
            
            meta.spike(s).adj = adj;
            meta.spike(s).windows = windows;
            meta.definitions.windows = {'non-spike','spikey','slow'};
            
            % Save the meta file after each spike run
            save(meta_file,'meta');
            
            t = toc;
            fprintf('Spike %d took %1.1f minutes.\n\n',s,t/60);
            
            if 0
                figure
                set(gcf,'position',[200 250 1175 400]);
                for tt = 1:size(windows,1)
                    subplot(1,size(windows,1),tt);
                    imagesc(squeeze(adj(2).adj(tt,:,:)));
                    colorbar
                    title(sprintf('%d s',tt))
                end
                
                beep
                pause
                close(gcf)
            end
            
        end
        
    end
    
end

end