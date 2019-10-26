function spike_networks(whichPts,do_simple_corr)

%% Parameters
merge = 1; % merge with existing?
do_car = 1;
pre_whiten = 0;
time_window = 0.5; %in seconds
n_chunks = 22; % number of time windows

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
    if do_simple_corr == 1
        out_folder = [pt_folder,'adj_simple/'];
    else
        out_folder = [pt_folder,'adj_coherence/'];
    end
    
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
            fprintf('Doing file %d of %d...\n',f,10);
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
            %old_values = values;
            values = pre_processing(values,do_car,pre_whiten);
            nchs = size(values,2);
            
            %% Figure out times for which I will be calculating adjacencies
            % The peak should be the very center of each file
            peak = round(size(values,1)/2);
            
            % Get index windows
            index_windows = zeros(n_chunks,2);
            tick_window = time_window*data.fs;
            
            if 0 % old method of defining windows
            
                % Get the time that I will call my spike window (that I will
                % remove when doing analysis)
                spike_window = peak + [-tick_window/2 tick_window/2];

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
                    continue
                end

                index_windows(ceil(n_chunks/2),:) = spike_window;

                for i = 1:ceil(n_chunks/2)-1
                    index_windows(i,1) = spike_window(1) - tick_window*(ceil(n_chunks/2)-i);
                    index_windows(i,2) = index_windows(i,1) + tick_window;
                end

                for i = ceil(n_chunks/2)+1:n_chunks
                    index_windows(i,1) = spike_window(2) + tick_window*(i-ceil(n_chunks/2)-1);
                    index_windows(i,2) = index_windows(i,1) + tick_window;
                end
            
            else % new way of defining windows; spike is in between the two middle ones
                for i = 1:n_chunks
                    index_windows(i,1) = peak - tick_window*n_chunks/2 + tick_window*(i-1);
                    index_windows(i,2) = peak - tick_window*n_chunks/2 + tick_window*(i);
                end
                
            end
            
            %% Get adjacency matrices
            %fprintf('Calculating functional networks...\n');

            if do_simple_corr == 0
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
                
            elseif do_simple_corr == 1
                adj = zeros(n_chunks,nchs,nchs);
                for tt = 1:n_chunks
                    % get appropriate points
                    temp_values = values(round(index_windows(tt,1)):round(index_windows(tt,2)),:); 
                    adj(tt,:,:) = get_simple_corr(temp_values);
                end
                
                if 1
                    figure
                    set(gcf,'position',[1 500 1400 297])
                    nplots = 5;
                    ha = tight_subplot(1,nplots);
                    
                    for i = 1:nplots
                        axes(ha(i))
                        imagesc(squeeze(adj(i+8,:,:)))
                    end
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