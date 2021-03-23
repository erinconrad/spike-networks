function spike_time_power

%% Parameters
max_time_off = 0.3;

%% Get locations
locations = spike_network_files;
main_folder = locations.main_folder;
data_folder = [main_folder,'data/'];
eeg_folder =  [main_folder,'results/eeg_data/'];
pt_file = [data_folder,'spike_structures/pt.mat'];
script_folder = locations.script_folder;
addpath(genpath(script_folder));

out_folder = [main_folder,'results/spike_time_power/'];
if exist(out_folder) == 0
    mkdir(out_folder);
end

pt = load(pt_file); % will create a structure called "pt"
pt = pt.pt;

%% Get manual spike times
sp = get_manual_times_from_excel(0);

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
    
    
    % Skip it if name is empty
    if isempty(pt(whichPt).name) == 1, continue; end
    name = pt(whichPt).name;
    fprintf('\nDoing %s\n',name);
    
    out_file = [name,'_spike_time_power.mat'];
    
    spike = load([eeg_folder,sprintf('%s_eeg.mat',name)]);
    spike = spike.spike;
    values = spike(1).data;
    nchs = size(values,2);
    fs = spike(1).fs;
    
    spike_powers = nan(length(spike),nchs);
    
    for s = 1:length(spike)
        if isempty(spike(s).time) == 1, continue; end
        values = spike(s).data;
        
        %% Pre-processing
        % Parameters 2 and 3 indicate whether to do CAR and pre-whiten,
        % respectively; 4 is whether to do notch filter
        %fprintf('Doing pre-processing...\n');
        %old_values = values;
        values = pre_processing(values,0,0,1,fs);
        
        % Bandpass filter the data to get the spikey component
        np = 6;
        fc = [5 70];%[2 70];
        fn = fs/2;
        [B,A] = butter(np,fc/fn);

        hp_data = zeros(size(values));
        for ich = 1:nchs
            hp_data(:,ich) = filtfilt(B,A,values(:,ich));
            %hp_data(:,ich) = filter(bpFilt,data(:,ich));
        end
        
        % restrict times to 0.3 s before to 0.3 s after the written time
        % this is to find the involved channels and peak IED channel
        narrow_indices = round(size(values,1)/2 - fs*max_time_off):...
            round(size(values,1)/2 + fs*max_time_off);
        narrow_data = hp_data(narrow_indices,:);
        
        % store deviations for all channels
        all_peak_powers = zeros(nchs,1);

        for ich = 1:nchs

            % Find the power
            bl = median(hp_data(:,ich));
            power = (narrow_data(:,ich)-bl).^2;

            % Find the peak deviation from the baseline
            [peak_dev,peak_index] = max(power);

           
            % all devs
            all_peak_powers(ich) = peak_dev;
        end

        % Find the biggest deviation channel
        [~,biggest_dev_ch] = max(all_peak_powers);
        
        if 0
            plot(values(:,biggest_dev_ch))
            pause
        end
        
        % Store the powers
        spike_powers(s,:) = all_peak_powers;
    end
    
    save([out_folder,out_file],'spike_powers');
    
end

end