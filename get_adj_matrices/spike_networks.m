function spike_networks(whichPts)

%% Parameters
merge = 1; % merge with existing?
do_car = 1;
pre_whiten = 1;
time_window = 1; %1 second time window
n_chunks = 11;

% The time surrounding the spike peak that I will call the spike window
spike_window_times = [-0.2 0.8];
% NEED TO ADJUST THIS!

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
 
            %% Pre-processing
            % Parameters 2 and 3 indicate whether to do CAR and pre-whiten,
            % respectively
            %fprintf('Doing pre-processing...\n');
            old_values = values;
            values = pre_processing(values,do_car,pre_whiten);
            nchs = size(values,2);
            
            %% Figure out times for which I will be calculating adjacencies
            % The peak should be the very center of each file
            peak = round(size(values,1)/2);
            
            % Get index windows
            index_windows = zeros(n_chunks,2); %6 is the spike window
            tick_window = time_window*data.fs;
            
            % Get the time that I will call my spike window (that I will
            % remove when doing analysis)
            spike_window = peak + spike_window_times*data.fs;
            
            % Show the spike window
            if 0
                sp_data = old_values(:,is_sp_ch);
                figure
                plot(sp_data);
                hold on
                plot(values(:,is_sp_ch))
                plot([spike_window(1) spike_window(1)],get(gca,'ylim'),'k--')
                plot([spike_window(2) spike_window(2)],get(gca,'ylim'),'k--')
                plot([peak peak],get(gca,'ylim'),'b--');
                beep
                pause
                close(gcf)
            end
            
            for i = 1:size(index_windows,1)
                
                % if we're at the spike window (6), then the index window
                % is the spike window. If i is 1, then it's 5 seconds
                % before. If i is 11, it's 5 seconds after
                index_windows(i,:) = spike_window + tick_window*(i-6);
            end
            
            %% Get adjacency matrices
            %fprintf('Calculating functional networks...\n');

            % Initialize adjacency matrix
            for ff = 1:size(freq_bands,1)
                adj(ff).name = freq_names{ff};
                adj(ff).adj = zeros(n_chunks,nchs,nchs);
            end
            
            for tt = 1:n_chunks
                
                % get appropriate points
                temp_values = values(round(index_windows(tt,1)):round(index_windows(tt,2)),:); 
                
                % Get adjacency matrices
                t_adj = get_adj_matrices(temp_values,data.fs,freq_bands);
                
                for ff = 1:size(freq_bands,1)
                    adj(ff).adj(tt,:,:) = t_adj(ff).adj;
                end
                
            end
            
            meta.spike(s).adj = adj;
            
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