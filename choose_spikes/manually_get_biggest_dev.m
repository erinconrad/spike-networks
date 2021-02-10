function manually_get_biggest_dev

%% Parameters
overwrite=0;
surround_time = 1;
do_notch = 1; % notch filter?
do_car = 1; % common average reference?
pre_whiten = 0; % remove the AR(1) component for a pre-whitening step?

%% Get file locations, load spike times and pt structure
locations = spike_network_files;
main_folder = locations.main_folder;
results_folder = [main_folder,'results/'];
data_folder = [main_folder,'data/'];
script_folder = locations.script_folder;
addpath(genpath(script_folder));
addpath(genpath(locations.BCT));
pt_file = [data_folder,'spike_structures/pt.mat'];
bct_folder = locations.BCT;
addpath(genpath(bct_folder));

eeg_folder = [results_folder,'eeg_data/'];
spike_rise_folder = [results_folder,'biggest_dev/'];


if exist(spike_rise_folder,'dir') == 0
    mkdir(spike_rise_folder);
end

listing = dir([eeg_folder,'*','_eeg.mat']);

all_names = {};
     
% Loop through eeg files (contains both spike and not a spike data)
for i = 1:length(listing)
    
    filename = listing(i).name;
    name_sp = split(filename,'_');
    name = name_sp{1};
    
    % skip not a spike files
    if contains(filename,'not') == 1, continue; end
    
    [a,b] = ismember(name,all_names);
    if a == 1
        pt_idx = b;
    else
        all_names = [all_names;name];
        pt_idx = length(all_names);
    end
    
    rise(pt_idx).name = name;
    
    fprintf('\nDoing %s...\n',name);
    
    % load eeg data
    spike = load([eeg_folder,filename]);
    spike = spike.spike;
    n_spikes = length(spike);
    values = spike(1).data;
    nchs = size(values,2);
    fs = spike(1).fs;
    
    % Initialize output data
    early_file = [spike_rise_folder,sprintf('%s_rise.mat',name)];
    
    % Load it if it exists to see how much we've already done
    if overwrite == 0
        if exist(early_file,'file') ~= 0
            early = load(early_file);
            early = early.early;

            % Find first unfinished spike
            start_spike = length(early.spike) + 1;
            
            fprintf('File already exists, loading and starting from unfinished spike.\n');

        else
            start_spike = 1;
            early.name = name;
        end
    else
        
        start_spike = 1;
        early.name = name;
    end
    
    % loop through spikes
    for s = start_spike:length(spike)
        
        fprintf('Doing spike %d of %d...\n',s,length(spike));
        
        % get eeg data
        data = spike(s).data; % ntimes x nch
        
        % pre process
        data = pre_processing(data,do_car,pre_whiten,do_notch,fs);
        
        % restrict to involved channels
        data_spike = data;
        
        % plot the spikes
        if 1
        figure
        set(gcf,'position',[-1343  194  500 804])
        offset = 0;
        ch_offsets = zeros(size(data_spike,2),1);
        ch_bl = zeros(size(data_spike,2),1);
        for ich = 1:size(data_spike,2)
            plot(linspace(-3,3,size(data_spike,1)),data_spike(:,ich)+offset);
            ch_offsets(ich) = offset;
            ch_bl(ich) = offset + median(data_spike(:,ich));
            hold on
            text(surround_time+0.05,ch_bl(ich),sprintf('%s',spike(s).chLabels{ich}))
            if ich<size(data_spike,2)
                offset = offset + max(data_spike(:,ich)) - min(data_spike(:,ich+1));
            end
        end
        
        xlim([-surround_time,surround_time]);
        
        title(sprintf('%s spike %d of %d',name,s,n_spikes),'fontsize',15)
        fprintf('Select the channel with the biggest spike deviation and press Enter\n');
        [~,y] = ginput;
        
        % Find closest channel
        [~,cl_ch] = (min(abs(ch_bl-y)));
        
        % Replot it to confirm
        hold off
        offset = 0;
        for ich = 1:size(data_spike,2)
            ch_offsets(ich) = offset;
            ch_bl(ich) = offset + median(data_spike(:,ich));
            if ich == cl_ch
                plot(linspace(-3,3,size(data_spike,1)),data_spike(:,ich)+offset,'linewidth',3);
                text(surround_time+0.05,ch_bl(ich),sprintf('%s',spike(s).chLabels{ich}))
            else
                plot(linspace(-3,3,size(data_spike,1)),data_spike(:,ich)+offset);
                text(surround_time+0.05,ch_bl(ich),spike(s).chLabels{ich})
            end

            xlim([-surround_time,surround_time]);
            if ich<size(data_spike,2)
                offset = offset + max(data_spike(:,ich)) - min(data_spike(:,ich+1));
            end
            hold on
        end

        inpt = input('\nIs this the correct channel? (y/n)\n','s');
        if strcmp(inpt,'y') || strcmp(inpt,'yes')
            biggest_dev_ch = cl_ch;
            fprintf('Designating %s\n',spike(s).chLabels{cl_ch});
        end
        
        early.spike(s).dev_ch = biggest_dev_ch;
        early.spike(s).dev_ch_label = spike(s).chLabels{cl_ch};
        
        % Save the meta file after each spike run
        save(early_file,'early');
        close(gcf)
        end
        
        
    end
    clear early
    
end