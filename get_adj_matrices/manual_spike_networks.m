function manual_spike_networks(whichPts,overwrite,do_simple_corr,time_window,not_a_spike)

%{
This function calculates functional networks for EEG data surrounding
manually detected spikes
%}

%{
Pt names
6 HUP074
7 HUP075
8 HUP078
9 HUP080
10 HUP082
11 HUP083
15 HUP094
16 HUP105
17 HUP106
18 HUP107
20 HUP116
23 Study017
24 Study019
30 Study028
31 Study029

%}

%% Parameters
do_notch = 1; % notch filter?
do_car = 1; % common average reference?
pre_whiten = 0; % remove the AR(1) component for a pre-whitening step?

freq_bands = [0.5 12;... %delta/theta/alpha
    12 25;... %beta
    30 40;... % low gamma
    96 106;... % high gamma
    106 256;... %ultra-high
    0.5 256;... %broadband    
    ]; 
freq_names = {'delta/theta/alpha','beta','low_gamma',...
    'high_gamma','ultra_high','broadband'};

%{
freq_bands = [0 4;... %delta
    4 8;...%theta
    8 12;...% alpha
    12 24;... %beta
    30 40;... % low gamma
    96 106;... % high gamma
    106 256;... %ultra-high
    0 256;... %broadband    
    ]; 


freq_names = {'delta','theta','alpha','beta','low_gamma',...
    'high_gamma','ultra_high','broadband'};
    %}
pos_text = '';
if length(time_window) == 1
    time_text = sprintf('%1.1f/',time_window);
else
    ntimes = length(time_window);
    all_times = time_window;
    true_window = all_times(2)-all_times(1);
    time_text = sprintf('%1.1f/',true_window);
    if time_window(1) == 0
        pos_text = 'pos_';
    end
end

if not_a_spike
    not_a_spike_text = '_not_spike_';
else
    not_a_spike_text = '_';
end

%% Get file locations, load spike times and pt structure
locations = spike_network_files;
main_folder = locations.main_folder;
eeg_folder = [main_folder,'results/eeg_data/'];
results_folder = [main_folder,'results/adj_mat/manual/'];
data_folder = [main_folder,'data/'];
script_folder = locations.script_folder;
addpath(genpath(script_folder));
pt_file = [data_folder,'spike_structures/pt.mat'];

pt = load(pt_file); % will create a structure called "pt"
pt = pt.pt;

%% Get manual spike times
sp = get_manual_times_from_excel(not_a_spike);

%{
% old spikes
sp = load([sp_folder,'sp.mat']);
sp = sp.sp;
%}

if isempty(whichPts) == 1
    whichPts = [];
    for i = 1:length(sp)
        if isempty(sp(i).name) == 0
            if sp(i).complete == 1
                whichPts = [whichPts,i];
            end
        end
    end
end

% Loop through patients
for whichPt = whichPts
    
    %% Prep patient
    
    % Skip it if name is empty
    if isempty(pt(whichPt).name) == 1, continue; end
    name = pt(whichPt).name;
    fprintf('\nDoing %s\n',name);
    
    % output folder
    if do_simple_corr == 1
        out_folder = [results_folder,'adj_simple/',pos_text,time_text];
    else
        out_folder = [results_folder,'adj_coherence/',pos_text,time_text];
    end
    
    if exist(out_folder,'dir') == 0
        mkdir(out_folder);
    end
 
    
    spike = load([eeg_folder,sprintf('%s%seeg.mat',name,not_a_spike_text)]);
    spike = spike.spike;

    % Initialize output data
    meta_file = [out_folder,sprintf('%s%sadj.mat',name,not_a_spike_text)];

    % Load it if it exists to see how much we've already done
    if overwrite == 0
        if exist(meta_file,'file') ~= 0
            meta = load(meta_file);
            meta = meta.meta;

            % Find first unfinished spike
            start_spike = length(meta.spike) + 1;
            
            fprintf('File already exists, loading and starting from unfinished spike.\n');
        end
    elseif overwrite == 1
        
        start_spike = 1;
        meta.name = name;
    elseif overwrite == 2
        meta = load(meta_file);
        meta = meta.meta;
        start_spike = 1;
    end
          
    % Loop through spikes
    for s = start_spike:length(spike)
        fprintf('Doing spike %d of %d...\n',s,length(spike));
        tic
        if isempty(spike(s).time) == 1, continue; end

        % Grab data
        values = spike(s).data;
        involved = spike(s).involved;
        chLabels = spike(s).chLabels;
        fs = spike(s).fs;
        meta.fs = fs;
        meta.chLabels = chLabels;
        meta.spike(s).time =spike(s).time;
        meta.spike(s).is_sp_ch = involved;
        
        %{
        if overwrite == 2
            if ~isequal(meta.freq_bands,freq_bands)
                error('Frequency bands do not line up.');
            end
        end
        meta.freq_bands = freq_bands;
        %}

        %% Pre-processing
        % Parameters 2 and 3 indicate whether to do CAR and pre-whiten,
        % respectively; 4 is whether to do notch filter
        %fprintf('Doing pre-processing...\n');
        %old_values = values;
        values = pre_processing(values,do_car,pre_whiten,do_notch,fs);
        nchs = size(values,2);
        
        if length(time_window) == 1
            n_chunks = round(size(values,1)/fs/time_window); % old way, same total time
        else
            n_chunks = ntimes;
        end
        %n_chunks = round(size(values,1)/fs); % alternate way, change total time so that it's the same number of chunks

        %% Figure out times for which I will be calculating adjacencies
        % The peak should be the very center of each file
        peak = round(size(values,1)/2);

        % Get index windows
        if length(time_window) == 1
            index_windows = zeros(n_chunks,2);
            tick_window = time_window*fs;

            for i = 1:n_chunks
                index_windows(i,1) = peak - tick_window*n_chunks/2 + tick_window*(i-1);
                index_windows(i,2) = peak - tick_window*n_chunks/2 + tick_window*(i);
            end
        else
            index_windows = zeros(n_chunks,2);
            
            for i = 1:n_chunks
                index_windows(i,1) = peak + round(time_window(i)*fs);
                index_windows(i,2) = peak + round(time_window(i)*fs) + round(true_window*fs);
            end
        end
        
        % Fix the first and the last to make sure they don't become
        % negative or beyond the total size
        index_windows(1,1) = max(index_windows(1,1),1);
        index_windows(end,2) = min(index_windows(end,2),size(values,1));
        
        % Try to line up with old time windows
        if overwrite == 2
            old_index_windows = meta.spike(s).index_windows;
            old_adj = meta.spike(s).adj;
            
            if isequal(old_index_windows(end,:),index_windows(1,:))
                append = 1; % I am adding on to the end
                index_windows(1,:) = []; % remove the first index window because it's the same
                final_index_windows = [old_index_windows;index_windows];
                n_chunks = n_chunks -1;
            elseif isequal(old_index_windows(1,:),index_windows(end,:))
                append = 2; % I am adding on before the beginning
                index_windows(end,:) = []; % remove last index window
                final_index_windows = [index_windows;old_index_windows];
                n_chunks = n_chunks -1;
            else
                error('I am not sure how to append this to the existing data');
            end
        else
            append = 0;
            final_index_windows = index_windows;
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
                t_adj = get_adj_matrices(temp_values,fs,freq_bands);

                for ff = 1:size(freq_bands,1)
                    
                    if append == 1
                        adj(ff).adj = [old_adj(ff).adj;zeros(n_chunks,nchs,nchs)];
                    elseif append == 2
                        adj(ff).adj = ...
                                [zeros(n_chunks,nchs,nchs);old_adj(ff).adj];
                    else
                        adj(ff).adj = zeros(n_chunks,nchs,nchs);
                    end
                    
                    if append == 1
                        % add onto old index windows
                        adj(ff).adj(tt+size(old_index_windows,1),:,:) = t_adj(ff).adj;
                    elseif append == 2
                        adj(ff).adj(tt,:,:) = t_adj(ff).adj;
                    else
                        adj(ff).adj(tt,:,:) = t_adj(ff).adj;
                    end
                end

            end

        elseif do_simple_corr == 1
            
            if append == 1
                adj(1).adj = [old_adj.adj;zeros(n_chunks,nchs,nchs)];
            elseif append == 2
                adj(1).adj = ...
                        [zeros(n_chunks,nchs,nchs);old_adj.adj];
            else
                adj(1).adj = zeros(n_chunks,nchs,nchs);
            end
            for tt = 1:n_chunks
                % get appropriate points
                temp_values = values(round(index_windows(tt,1)):round(index_windows(tt,2)),:); 
                if append == 1
                    adj(1).adj(tt+size(old_index_windows,1),:,:) = get_simple_corr(temp_values);
                elseif append == 2
                    adj(1).adj(tt,:,:) = get_simple_corr(temp_values);
                else
                    adj(1).adj(tt,:,:) = get_simple_corr(temp_values);
                end
            end

            if 0
                figure
                set(gcf,'position',[1 500 1400 297])
                nplots = 5;
                ha = tight_subplot(1,nplots);

                for i = 1:nplots
                    axes(ha(i))
                    imagesc(squeeze(adj(1).adj(i+8,:,:)))
                end
            end
        end

        meta.spike(s).adj = adj;
        meta.spike(s).index_windows = final_index_windows;
        
        error('look');

        % Save the meta file after each spike run
        save(meta_file,'meta');

        t = toc;
        fprintf('Spike %d took %1.1f minutes.\n\n',s,t/60);

        if 0
            figure

            set(gcf,'position',[200 250 1175 400]);
            for tt = 1:n_chunks-1
                subplot(2,5,tt);
                imagesc(squeeze(adj(8).adj(tt,:,:)));
                colorbar
                title(sprintf('%d s',tt))
            end
            beep
            %return
            pause
            close(gcf)
        end
    end

    % Once I've done all spikes for the file, make a finished.mat
   % save([out_folder,sprintf('finished.mat')],'finished');
        
    clear meta
end


end