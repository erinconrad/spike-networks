function get_pdc(whichPts,overwrite,time_window,not_a_spike)


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

ntimes = length(time_window);
all_times = time_window;
true_window = all_times(2)-all_times(1);
time_text = sprintf('%1.1f/',true_window);

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

for whichPt = whichPts
    
    %% Prep patient
    
    % Skip it if name is empty
    if isempty(pt(whichPt).name) == 1, continue; end
    name = pt(whichPt).name;
    fprintf('\nDoing %s\n',name);

    out_folder = [results_folder,'pdc/',time_text];

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
        else
            start_spike = 1;
            clear meta
            meta.name = name;
        end
    elseif overwrite == 1
        
        start_spike = 1;
        meta.name = name;
    elseif overwrite == 2
        meta = load(meta_file);
        meta = meta.meta;
        start_spike = 1;
    end

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
        meta.spike(s).biggest_dev = spike(s).biggest_dev;


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

        index_windows = zeros(n_chunks,2);
            
        for i = 1:n_chunks
            index_windows(i,1) = peak + round(time_window(i)*fs);
            index_windows(i,2) = peak + round(time_window(i)*fs) + round(true_window*fs);
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
                fprintf('Warning, I am not sure how to append this to the existing data. Skipping spike\n');
                continue
            end
        else
            append = 0;
            final_index_windows = index_windows;
        end

        pdc_out = pdc(values,fs,index_windows,freq_bands);
        for tt = 1:n_chunks
            for ff = 1:size(freq_bands,1)
                adj(ff).adj(tt,:,:) = squeeze(pdc_out(:,:,ff,tt));
            end
        end

        meta.spike(s).adj = adj;
        meta.spike(s).index_windows = final_index_windows;

        % Save the meta file after each spike run
        save(meta_file,'meta');

        t = toc;
        fprintf('Spike %d took %1.1f minutes.\n\n',s,t/60);

    end


    clear meta
end


end