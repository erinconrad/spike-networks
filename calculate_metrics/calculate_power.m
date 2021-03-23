function calculate_power(time_window,not_a_spike)

%% Description
%{
This function calculates frequency specific power
%}

%% Parameters
do_notch = 1; % notch filter?
do_car = 1; % common average reference?
pre_whiten = 0; % remove the AR(1) component for a pre-whitening step?


%% Proposed frequencies, following Ren et al., Neurology, 2015
freq_bands = [0.5 30;... %sub-gamma
    30 100;... % ;low gamma
    100 256]; % high_gamma

freq_names = {'sub-gamma','low gamma','high gamma'};

n_f = size(freq_bands,1);

%% Time windows
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
results_folder = [main_folder,'results/'];
data_folder = [main_folder,'data/'];
script_folder = locations.script_folder;
addpath(genpath(script_folder));
pt_file = [data_folder,'spike_structures/pt.mat'];
output_folder = [results_folder,'all_powers/',time_text,'/'];

if exist(output_folder,'dir') == 0
    mkdir(output_folder);
end

pt = load(pt_file); % will create a structure called "pt"
pt = pt.pt;

%% Get manual spike times
sp = get_manual_times_from_excel(not_a_spike);

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
    clear power
    
    % Skip it if name is empty
    if isempty(pt(whichPt).name) == 1, continue; end
    name = pt(whichPt).name;
    fprintf('\nDoing %s\n',name);
    
    spike = load([eeg_folder,sprintf('%s%seeg.mat',name,not_a_spike_text)]);
    spike = spike.spike;
    n_spikes = length(spike);
    values = spike(1).data;
    nchs = size(values,2);
    fs = spike(1).fs;
    
 
    n_windows = ntimes;

    
    out_file = [output_folder,sprintf('%s%spower.mat',name,not_a_spike_text)];
    
    % Initialize arrays
    abs_power_array = nan(n_spikes,n_windows,nchs);
    ers_array = nan(n_spikes,n_windows,nchs,n_f);
    
    for s = 1:length(spike)
        if isempty(spike(s).time) == 1, continue; end
        values = spike(s).data;
        
        %% Pre-processing
        % Parameters 2 and 3 indicate whether to do CAR and pre-whiten,
        % respectively; 4 is whether to do notch filter
        %fprintf('Doing pre-processing...\n');
        %old_values = values;
        values = pre_processing(values,do_car,pre_whiten,do_notch,fs);
        
        % subtract baseline
        X = values-median(values,1);
        
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
        index_windows = zeros(n_windows,2);

        for i = 1:n_windows
            index_windows(i,1) = peak + round(time_window(i)*fs);
            index_windows(i,2) = peak + round(time_window(i)*fs) + round(true_window*fs);
        end
        
        % Fix the first and the last to make sure they don't become
        % negative or beyond the total size
        index_windows(1,1) = max(index_windows(1,1),1);
        index_windows(end,2) = min(index_windows(end,2),size(values,1));
        
        %% Get ERS
        
        
        % Get ERS
        for ich = 1:size(X,2)
            for t = 1:n_windows
                
                % Restrict to specific channel and time window
                Xtemp = X(max(1,round(index_windows(t,1))):...
                    min(size(X,1),round(index_windows(t,2))),ich);
                
                % Get the frequency-specific powers
                powers = get_power(Xtemp,fs,freq_bands);
                ers_array(s,t,ich,:) = powers;
            end
        end
        
        
        
        
        %% Get absolute power
        
        dev_squared = X.^2;
        
        % Loop through chs
        for ich = 1:size(X,2)
            dev_squared_ch = dev_squared(:,ich);
            
            % now, get the average power in each time window for that spike
            for t = 1:n_windows
                abs_power_array(s,t,ich) = mean(dev_squared_ch(max(1,round(index_windows(t,1)))...
                    :min(length(dev_squared_ch),round(index_windows(t,2)))));

            end
        end
        
        
    end
    
    %% Fill struct
    power.ers = ers_array;
    power.nchs = nchs;
    power.abs_power = abs_power_array;
    power.name = name;
    power.time_window = time_window;
    power.index_windows = index_windows;
    power.n_windows = n_windows;
    power.freq_names = freq_names;
    power.freq_bands = freq_bands;
    power.fs = fs;
    
    % save
    save(out_file,'power');
end

end