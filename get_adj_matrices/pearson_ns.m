function pearson_ns(whichPts,time_window,not_a_spike)

%% Parameters
do_notch = 1; % notch filter?
do_car = 1; % common average reference?
pre_whiten = 0; % remove the AR(1) component for a pre-whitening step?

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
bct_folder = locations.BCT;
addpath(genpath(bct_folder));

pt_file = [data_folder,'spike_structures/pt.mat'];

pt = load(pt_file); % will create a structure called "pt"
pt = pt.pt;

%% Get manual spike times
sp = get_manual_times_from_excel(not_a_spike);


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
    
    out_folder = [results_folder,'simple_ns/',time_text];
    
    if exist(out_folder,'dir') == 0
        mkdir(out_folder);
    end
 
    
    spike = load([eeg_folder,sprintf('%s%seeg.mat',name,not_a_spike_text)]);
    spike = spike.spike;
    
    % Initialize output data
    out_file = [out_folder,sprintf('%s%sadj.mat',name,not_a_spike_text)];

    
    clear meta
    
    % get basic info
    fs = spike(1).fs;
    nchs = length(spike(1).chLabels);
    meta.name = name;
    meta.fs = fs;
    meta.nchs = nchs;
    meta.time_window = time_window;
    n_spikes = length(spike);
    n_times = length(time_window);
    
    ns = nan(n_spikes,n_times,nchs);
        
    % Loop through spikes
    for s = 1:length(spike)
        fprintf('Doing spike %d of %d...\n',s,length(spike));
        tic
        if isempty(spike(s).time) == 1, continue; end

        % Grab data
        values = spike(s).data;
        fs = spike(s).fs;
        
       
        %% Pre-processing
        % Parameters 2 and 3 indicate whether to do CAR and pre-whiten,
        % respectively; 4 is whether to do notch filter
        %fprintf('Doing pre-processing...\n');
        %old_values = values;
        values = pre_processing(values,do_car,pre_whiten,do_notch,fs);
        
        n_chunks = ntimes;
        peak = round(size(values,1)/2);
        index_windows = zeros(n_chunks,2);
            
        for i = 1:n_chunks
            index_windows(i,1) = peak + round(time_window(i)*fs);
            index_windows(i,2) = peak + round(time_window(i)*fs) + round(true_window*fs);
        end
        index_windows(1,1) = max(index_windows(1,1),1);
        index_windows(end,2) = min(index_windows(end,2),size(values,1));
        meta.index_windows = index_windows;
        
        
        for tt = 1:n_chunks
            % get appropriate points
            temp_values = values(round(index_windows(tt,1)):round(index_windows(tt,2)),:); 
            
            % Calculate adjacency matrix based on Pearson correlation
            % coefficient
            adj = get_simple_corr(temp_values);
            
            % Calculate node strength
            ns_temp = strengths_und(adj);
            ns(s,tt,:) = ns_temp;
        end
        
    end
    meta.data = ns;
    % save
    save(out_file,'meta');
    
end

end