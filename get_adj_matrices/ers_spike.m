function ers_spike(overwrite,time_window)

%% Parameters
do_notch = 1; % notch filter?
do_car = 1; % common average reference?
pre_whiten = 0; % remove the AR(1) component for a pre-whitening step?

freq_bands = [0.5 4;... %delta
    4 8;...%theta
    8 12;...% alpha
    12 24;... %beta
    30 40;... % low gamma
    96 106;... % high gamma
    106 256;... %ultra-high
    0.5 256;... %broadband    
    ]; 
freq_names = {'delta','theta','alpha','beta','low_gamma',...
    'high_gamma','ultra_high','broadband'};
n_f = size(freq_bands,1);
time_text = sprintf('%1.1f/',time_window);

%% Get file locations, load spike times and pt structure
locations = spike_network_files;
main_folder = locations.main_folder;
eeg_folder = [main_folder,'results/eeg_data/'];
results_folder = [main_folder,'results/'];
data_folder = [main_folder,'data/'];
script_folder = locations.script_folder;
addpath(genpath(script_folder));
pt_file = [data_folder,'spike_structures/pt.mat'];
output_folder = [results_folder,'ers/',time_text,'/'];
if exist(output_folder,'dir') == 0
    mkdir(output_folder);
end

pt = load(pt_file); % will create a structure called "pt"
pt = pt.pt;

%% Get manual spike times
sp = get_manual_times_from_excel;

whichPts = [];
for i = 1:length(sp)
    if isempty(sp(i).name) == 0
        if sp(i).complete == 1
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
    
    spike = load([eeg_folder,sprintf('%s_eeg.mat',name)]);
    spike = spike.spike;
    n_spikes = length(spike);
    values = spike(1).data;
    nchs = size(values,2);
    fs = spike(1).fs;
    n_windows = round(size(values,1)/fs/time_window);
    
    out_file = [output_folder,sprintf('%s_ers.mat',name)];
    % Load it if it exists to see how much we've already done
    if overwrite == 0
        if exist(out_file,'file') ~= 0
            fprintf('Already done %s, skipping...\n',name);
            continue;
        end
    end
    
    % Initialize array of ERS
    ers_array = nan(n_spikes,n_windows,n_f,nchs);
    ers_involved = nan(n_spikes,n_windows,n_f);
    
    for s = 1:length(spike)
        if isempty(spike(s).time) == 1, continue; end
        values = spike(s).data;
        involved = spike(s).involved;
        
        %% Pre-processing
        % Parameters 2 and 3 indicate whether to do CAR and pre-whiten,
        % respectively; 4 is whether to do notch filter
        %fprintf('Doing pre-processing...\n');
        %old_values = values;
        values = pre_processing(values,do_car,pre_whiten,do_notch,fs);
        
        
        % show spike
        if 0
        involved_chs = find(involved);
        figure
        set(gcf,'position',[150 450 1000 350])
        for ich = 1:length(involved_chs)
            ch = involved_chs(ich);
            plot(linspace(-spike(1).surround_time,spike(1).surround_time,...
                size(values,1)),values(:,ch))
            title(sprintf('Ch %d',ch))
            pause
        end
        end
        

        
        if 0
        % Plot the spectrogram
        x = values(:,77);
        spectrogram(x,128,120,128,fs,'yaxis')
        pause
        close(gcf)
        end
        
         %% Figure out times for which I will be calculating ERS
        % The peak should be the very center of each file
        peak = round(size(values,1)/2);

        % Get index windows
        index_windows = zeros(n_windows,2);
        tick_window = time_window*fs;
        
        for i = 1:n_windows
            index_windows(i,1) = peak - tick_window*n_windows/2 + tick_window*(i-1);
            index_windows(i,2) = peak - tick_window*n_windows/2 + tick_window*(i);
        end
        
        % Fix the first and the last to make sure they don't become
        % negative or beyond the total size
        index_windows(1,1) = max(index_windows(1,1),1);
        index_windows(end,2) = min(index_windows(end,2),size(values,1));
        
        % Get ERS
        for ich = 1:nchs
            for t = 1:n_windows
                X = values(round(index_windows(t,1)):round(index_windows(t,2)),ich);
                powers = get_power(X,fs,freq_bands);
                ers_array(s,t,:,ich) = powers;
            end
        end
        
        ers.spike(s).involved = involved;
        ers.spikes(s).index_windows = index_windows;
        ers.spike(s).ers_involved = nanmean(ers_array(s,:,:,involved),4);
        ers_involved(s,:,:) = nanmean(ers_array(s,:,:,involved),4);
        
    end
    
    % Fill struct
    ers.name = name;
    ers.time_window = time_window;
    ers.n_windows = n_windows;
    ers.powers = ers_array;
    ers.powers_avg_involved = squeeze(nanmean(ers_involved,1));
    ers.powers_involved = ers_involved;
    ers.freq_names = freq_names;
    ers.freq_bands = freq_bands;
    
    % save
    save(out_file,'ers');
end

end