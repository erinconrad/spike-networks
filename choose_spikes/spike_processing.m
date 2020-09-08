function spike_processing(whichPts,do_save,overwrite,not_a_spike)

%{
This function stores eeg data for manually detected spikes
%}

%% Parameters
surround_time = 3; % how many seconds before and after each spike to take for initial analysis
thresh_spike = 8; % How high above baseline the filtered data must be

%% Get locations
locations = spike_network_files;
script_folder = locations.script_folder;
main_folder = locations.main_folder;
pwname = locations.pwfile;
addpath(genpath(script_folder));

if isempty(locations.ieeg_folder) == 0
    addpath(genpath(locations.ieeg_folder));
end

results_folder =  [main_folder,'results/eeg_data/'];
if exist(results_folder,'dir') == 0
    mkdir(results_folder);
end

pt_folder = [main_folder,'data/spike_structures/'];
pt = load([pt_folder,'pt.mat']);
pt = pt.pt;

if not_a_spike
    not_a_spike_text = '_not_spike_';
else
    not_a_spike_text = '_';
end

%% Get the spike times
sp_folder = [main_folder,'data/manual_spikes/'];
sp = get_manual_times_from_excel(not_a_spike);

% Get the patients
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
    
    clear spike
    
    name = sp(whichPt).name;
    
    if isempty(name) == 1
        continue;
    elseif sp(whichPt).complete == 0 % Skip if we don't have complete spike data
        continue;
    else
        fprintf('Doing %s...\n',name);
    end
    
    % Skip if I already did it
    if overwrite == 0
        if exist([results_folder,name,not_a_spike_text,'eeg.mat'],'file') ~= 0
            fprintf('Already did %s, skipping...\n',name);
            continue;
        end
    end
    
    n_spikes = length(sp(whichPt).spike);
    
    % get ieeg name and sampling rate
    ieeg_name = pt(whichPt).ieeg_name;
    fs = pt(whichPt).fs;
    
    % get which channels to do: this tells me which ieeg channels to take
    % and in what order so that they should perfectly align with the
    % electrode designations in pt.new_elecs
    chs = pt(whichPt).new_elecs.ch_order;
    
    % Download dummy data
    data = download_eeg(ieeg_name,[],pwname,1,[]);
    fs = data.fs;
    chLabels = data.chLabels(:,1);
    chLabels = cellfun(@ieeg_ch_parser,chLabels,'UniformOutput',false);
    
    % Compare chLabels to pt(whichPt).new_elecs.names (should be identical
    % after parsing)
    if isequal(chLabels(chs),pt(whichPt).new_elecs.names) == 0
        error('what\n');
    end
    
    for s = 1:n_spikes
        
        %% Load the data surrounding the spike
        s_time = sp(whichPt).spike(s);
        which_times = [s_time-surround_time,s_time+surround_time];
        
        % get indices
        indices = round(which_times(1)*fs):round(which_times(2)*fs);

        % Get the data
        data = download_eeg(ieeg_name,indices,pwname,0,chs); 
        data = data.values;
        nch = size(data,2);
        
        % Store the data
        spike(s).data = data;
        spike(s).chLabels = chLabels(chs);
        spike(s).fs = fs;
        spike(s).times = which_times;
        spike(s).surround_time = surround_time;
        spike(s).thresh_spike = thresh_spike;
        
        %% Find peak spike time for each channel and which channels I will say are involved in the spike
        
        % alternate bandpass filter
        %{
        bpFilt = designfilt('bandpassiir','FilterOrder',6, ...
         'HalfPowerFrequency1',2,'HalfPowerFrequency2',70, ...
         'SampleRate',fs);
        %}
        
            % Common average reference and notch filter the data
            data = pre_processing(data,1,0,1,fs);

            % Bandpass filter the data to get the spikey component
            np = 6;
            fc = [5 70];%[2 70];
            fn = fs/2;
            [B,A] = butter(np,fc/fn);

            hp_data = zeros(size(data));
            for ich = 1:nch
                hp_data(:,ich) = filtfilt(B,A,data(:,ich));
                %hp_data(:,ich) = filter(bpFilt,data(:,ich));
            end

            if 0
                figure
                for tch = 1:nch
                    plot(data(:,tch));
                    hold on
                    plot(hp_data(:,tch));
                    title(sprintf('%d %s',tch,chLabels{tch}))
                    pause
                    close(gcf)
                end
            end


            % restrict times to 1 s before to 1 s after the written time
            narrow_indices = max(round(length(indices)/2 - fs),1):...
                min(round(length(indices)/2 + fs),length(indices));
            narrow_data = hp_data(narrow_indices,:);
            narrow_unfiltered = data(narrow_indices,:);


            % initialize additional data for storage
            spike(s).peak_time = zeros(nch,1);
            spike(s).rel_dev = zeros(nch,1);
            spike(s).involved = zeros(nch,1);

            % store deviations for all channels
            all_devs = zeros(nch,1);

            for ich = 1:nch

                % Find the deviation from the baseline
                bl = median(narrow_data(:,ich));
                dev = abs(narrow_data(:,ich)-bl);

                % Find the peak deviation from the baseline
                [peak_dev,peak_index] = max(dev);

                % Find the peak time (adjust for index)
                peak_time = (peak_index + narrow_indices(1) -1 + indices(1) -1)/fs;

                % Get the relative deviation in relation to standard deviation
                % of the deviation
                rel_dev = peak_dev/std(dev);

                % Decide if the channel is involved in the spike
                involved = rel_dev > thresh_spike;

                % all devs
                all_devs(ich) = rel_dev;

                % Store this
                spike(s).rel_dev(ich) = rel_dev;
                spike(s).peak_time(ich) = peak_time;
                spike(s).involved(ich) = involved;

                if 0
                    plot(narrow_unfiltered(:,ich))
                    hold on
                    plot(narrow_data(:,ich))
                    xl = get(gca,'xlim');
                    plot([xl(1) xl(2)],[median(narrow_data(:,ich)) median(narrow_data(:,ich))])
                    if involved == 1
                        plot(peak_index,narrow_data(peak_index,ich),'bo','MarkerSize',30)
                    else
                        plot(peak_index,narrow_data(peak_index,ich),'rx','MarkerSize',30)
                    end
                    title(chLabels(ich))
                    pause
                    close(gcf)
                end

            end

            %% Order the involved channels by latency
            spike(s).involved = logical(spike(s).involved);
            in_chs = find(spike(s).involved);
            [sorted_times,I] = sort(spike(s).peak_time(spike(s).involved));
            spike(s).ordered_chs = [in_chs(I),sorted_times];

            %% Find the biggest deviation channel
            [~,biggest_dev_ch] = max(all_devs);
            spike(s).biggest_dev = biggest_dev_ch;


            % Plot
            if 0
                figure
                set(gcf,'position',[100 100 1000 500])
                offset = 0;
                for ich = 1:size(spike(s).ordered_chs,1)
                    plot(data(:,spike(s).ordered_chs(ich,1))-offset)
                    hold on
                    offset = offset - 400;
                end
                pause
                close(gcf)
            end

            if 0
                figure
                set(gcf,'position',[100 100 1000 500])
                plot(data(:,spike(s).biggest_dev))
                pause
                close(gcf)
            end

            %% Set the spike "time" to be the median of the peak times
            spike(s).time = median(sorted_times);
        
        
    end
    if do_save == 1
        
        save([results_folder,name,not_a_spike_text,'eeg.mat'],'spike');
        
            
    end
    
end


end