function spike_networks(whichPts)

%% Parameters
do_car = 1;
pre_whiten = 0;
time_window = 1; %1 second time window

% The time prior to the spike peak that will start the 
pre_spike_time = 0.05;
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

times = load(spike_times_file); % will result in a structure called "times"
times = times.times;
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
        fprintf('Doing file %d of %d...\n',f,10);
        spike = load([pt_folder,sprintf('spikes_%d.mat',f)]);
        spike = spike.spike;
        
        % Loop through spikes
        for s = 1:length(spike)
            fprintf('Doing spike %d of %d...\n',s,length(spike));
            tic
            if isempty(spike(s).time) == 1, continue; end
            
            
            
            % Grab the appropriate channels
            values = spike(s).values(:,~ignore);
            
            % test plot of spike
            if 0
                
                % Get spike ch
                is_sp_ch = strcmp(ch_labels,spike(s).label);
                sp_data = spike(s).values(:,is_sp_ch);
                figure
                plot(sp_data);
                pause
                close(gcf)
            end
            
            %% Pre-processing
            % Parameters 2 and 3 indicate whether to do CAR and pre-whiten,
            % respectively
            %fprintf('Doing pre-processing...\n');
            values = pre_processing(values,do_car,pre_whiten);
            nchs = size(values,2);
            
            %% Get adjacency matrices
            %fprintf('Calculating functional networks...\n');
            % Break the data up into chunks
            tick_window = time_window*data.fs;
            n_chunks = floor(size(values,1)/tick_window);
            
            % Initialize adjacency matrix
            for ff = 1:size(freq_bands,1)
                adj(ff).name = freq_names{ff};
                adj(ff).adj = zeros(n_chunks,nchs,nchs);
            end
            
            for tt = 1:n_chunks
                
                % get appropriate points
                temp_values = values((tt-1)*tick_window+1:...
                    min(tt*tick_window,size(values,1)),:); % this will cut off last point for some, which is fine
                
                % Get adjacency matrices
                t_adj = get_adj_matrices(temp_values,data.fs,freq_bands);
                
                for ff = 1:size(freq_bands,1)
                    adj(ff).adj(tt,:,:) = t_adj(ff).adj;
                end
                
            end

            t = toc;
            fprintf('Spike %d took %1.1f minutes.\n\n',s,t/60);
            
            if 1
                figure
                set(gcf,'position',[200 250 1175 400]);
                for tt = 1:n_chunks
                    subplot(2,5,tt);
                    imagesc(squeeze(adj(4).adj(tt,:,:)));
                    colorbar
                    title(sprintf('%d s',tt))
                end
                pause
                close(gcf)
            end
        end
    end
    
end


end