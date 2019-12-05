function spike_processing(whichPts)

%% Parameters
surround_time = 6;
thresh_spike = 6;

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

for whichPt = whichPts
    
    name = sp(whichPt).name;
    
    if isempty(name) == 1
        continue;
    else
        fprintf('Doing %s...\n',name);
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
        s_time = sp(whichPt).spike(s).time;
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
        
        %% Find peak spike time for each channel and which channels I will say are involved in the spike
        
        % alternate bandpass filter
        %{
        bpFilt = designfilt('bandpassiir','FilterOrder',6, ...
         'HalfPowerFrequency1',2,'HalfPowerFrequency2',70, ...
         'SampleRate',fs);
        %}
        
        % Bandpass filter the data to get the spikey component
        np = 6;
        fc = [2 70];
        fn = fs/2;
        [B,A] = butter(np,fc/fn);
        
        hp_data = zeros(size(data));
        for ich = 1:nch
            hp_data(:,ich) = filtfilt(B,A,data(:,ich));
            %hp_data(:,ich) = filter(bpFilt,data(:,ich));
        end
        
        if 1
            figure
            tch = 20;
            plot(data(:,tch));
            hold on
            plot(hp_data(:,tch));
            pause
            close(gcf)
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

            
            % Store this
            spike(s).rel_dev(ich) = rel_dev;
            spike(s).peak_time(ich) = peak_time;
            spike(s).involved(ich) = involved;
            
            if 1
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
                pause
                close(gcf)
            end
            
        end
        
        %% Order the involved channels by latency
        spike(s).involved = logical(spike(s).involved);
        in_chs = find(spike(s).involved);
        [sorted_times,I] = sort(spike(s).peak_time(spike(s).involved));
        spike(s).ordered_chs = [in_chs(I),sorted_times];
        
        % Plot
        if 0
            figure
            set(gcf,'position',[100 100 1000 500])
            offset = 0;
            for ich = 1:length(spike(s).ordered_chs)
                plot(data(:,ich)-offset)
                hold on
                offset = offset - 400;
            end
            pause
            close(gcf)
        end
        
        %% Set the spike "time" to be the median of the peak times
        spike(s).time = median(sorted_times);
        
    end

    save([results_folder,name,'_eeg.mat'],'spike');
    
end


end