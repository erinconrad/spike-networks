function pdc_time(whichPts,overwrite,not_a_spike,which_times,which_freq)

%{
Consider not doing common average reference
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

    out_folder = [results_folder,'pdc_time/'];

    if exist(out_folder,'dir') == 0
        mkdir(out_folder);
    end
 
    
    spike = load([eeg_folder,sprintf('%s%seeg.mat',name,not_a_spike_text)]);
    spike = spike.spike;

    % Initialize output data
    meta_file = [out_folder,sprintf('%s%spdc_time.mat',name,not_a_spike_text)];

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
        dur = round(size(values,1)/fs * 1e2)/(1e2);
        mid = dur/2;
        
        % Restrict times to the desired times
        adj_times = which_times + mid;
        indices = round(adj_times(1)*fs:adj_times(2)*fs);
        
        curr_values = values(indices,:);
        
        [pdc_out,dtf_out] = pdc(curr_values,fs,[1 size(curr_values,1)],freq_bands(which_freq,:));
        meta.spike(s).pdc = pdc_out;
        meta.spike(s).dtf = dtf_out;
        meta.which_times = which_times;
        meta.spike(s).index_windows = indices;
        meta.which_freq = which_freq;
        meta.freq_name = freq_names{which_freq};
        
        % Save the meta file after each spike run
        save(meta_file,'meta');

        t = toc;
        fprintf('Spike %d took %1.1f minutes.\n\n',s,t/60);
        
    end
    clear meta

end

end