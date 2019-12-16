function manual_spike_networks(whichPts,do_simple_corr)

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
eeg_folder = [main_folder,'results/eeg_data/'];
results_folder = [main_folder,'results/adj_mat/'];
data_folder = [main_folder,'data/'];
script_folder = locations.script_folder;
addpath(genpath(script_folder));
pt_file = [data_folder,'spike_structures/pt.mat'];

pt = load(pt_file); % will create a structure called "pt"
pt = pt.pt;

sp_folder = [main_folder,'data/manual_spikes/'];
sp = load([sp_folder,'sp.mat']);
sp = sp.sp;

if isempty(whichPts) == 1
    whichPts = [];
    for i = 1:length(sp)
        if isempty(sp(i).name) == 0
            whichPts = [whichPts,i];
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
        out_folder = [results_folder,'adj_simple/'];
    else
        out_folder = [results_folder,'adj_coherence/'];
    end
    
    if exist(out_folder,'dir') == 0
        mkdir(out_folder);
    end

    
    spike = load([eeg_folder,sprintf('%s_eeg.mat',name)]);
    spike = spike.spike;

    % Initialize output data
    meta_file = [out_folder,sprintf('%s_adj.mat',name)];

    % Load it if it exists to see how much we've already done
    if exist(meta_file,'file') ~= 0
        meta = load(meta_file);
        meta = meta.meta;

        % Find first unfinished spike
        start_spike = length(meta.spike) + 1;

    else
        start_spike = 1;
        meta.name = name;
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
        time = spike(s).time;
        meta.fs = fs;
        meta.chLabels = chLabels;
        meta.spike(s).time =spike(s).time;
        meta.spike(s).is_sp_ch = involved;

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
        tick_window = time_window*fs;
        
        for i = 1:n_chunks
            index_windows(i,1) = peak - tick_window*n_chunks/2 + tick_window*(i-1);
            index_windows(i,2) = peak - tick_window*n_chunks/2 + tick_window*(i);
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
            adj(1).adj = zeros(n_chunks,nchs,nchs);
            for tt = 1:n_chunks
                % get appropriate points
                temp_values = values(round(index_windows(tt,1)):round(index_windows(tt,2)),:); 
                adj(1).adj(tt,:,:) = get_simple_corr(temp_values);
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
            %return
            pause
            close(gcf)
        end
    end

    % Once I've done all spikes for the file, make a finished.mat
   % save([out_folder,sprintf('finished.mat')],'finished');
        
    
end


end