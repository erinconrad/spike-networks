function out = get_metric_data(met,time_window)

%% Get file locations, load spike times and pt structure
locations = spike_network_files;
main_folder = locations.main_folder;
results_folder = [main_folder,'results/'];
power_folder = [results_folder,'all_powers/'];
ns_folder = [results_folder,'ns/'];

if contains(met,'ers') || contains(met,'power')
    main_folder = power_folder;
elseif contains(met,'ns')
    main_folder = ns_folder;
end

time_name = sprintf('%1.1f/',time_window);
time_folder = [main_folder,time_name];
pt_listing = dir([time_folder,'*.mat']);

%% Get patient names
all_names = {};
for i = 1:length(pt_listing)
    fname = pt_listing(i).name;
    pt_name = strsplit(fname,'_');
    pt_name = pt_name{1};
    if ~ismember(pt_name,all_names)
        all_names = [all_names;pt_name];
    end
end

% Loop through patients
for p = 1:length(all_names)
    name = all_names{p};
    
    %initialize spike and not spike listings
    sp_or_not_listing = cell(2,1);
    
    %% Find spike and not spike listings
    pt_listing = dir([time_folder,'*',name,'*']);
    if length(pt_listing) ~= 2, error('what'); end
    
    if contains(pt_listing(1).name,'not')
        % Assign the first one to the not a spike listing (the 2nd)
        sp_or_not_listing{2} = pt_listing(1).name;
        sp_or_not_listing{1} = pt_listing(2).name;
    else
        % Assign the first one to the spike listing (the first)
        sp_or_not_listing{1} = pt_listing(1).name;
        sp_or_not_listing{2} = pt_listing(2).name;
    end
    
    
    % Loop through spike and not spike
    for sp = 1:2
    
        %% Load the appropriate file
        info = load([time_folder,sp_or_not_listing{sp}]);
        
        %% Load up metric info
        switch met

            %% Get power
            case 'abs_power'
                index_windows = info.power.index_windows;
                data = info.power.abs_power;
                times = info.power.time_window;

            %% Get ERS
            case 'ers'
                index_windows = info.power.index_windows;
                data = info.power.ers;
                times = info.power.time_window;

            %% Get NS
            case 'ns'
                index_windows = info.metrics.index_windows;
                data = info.metrics.data;
                
                % Fix to get times (since I never stored them)
                fs = info.metrics.fs;
                index_diff = index_windows(1,2) - index_windows(1,1); 
                time_diff = round(index_diff*1e1/fs)/(1e1);
                first_time = round(index_windows(1,1)*1e1/fs)/(1e1);
                last_time = round(index_windows(end,1)*1e1/fs)/(1e1);
                old_times = first_time:time_diff:last_time;
                mid_time = 3; % middle of file is 3 seconds in
                times = old_times-mid_time;
                
        end
        
        %% Re-structure
        out.index_windows = index_windows;
        out.times = times;
        
        % Get array size
        sz = size(data);
        if length(sz) == 3
            nfreq = 1;
        else
            nfreq = 3;
        end
        
        for f = 1:nfreq
            out.freq(f).pt(p).name = name;
            
            if length(sz) == 3
                out.freq(f).pt(p).sp_or_not(sp).data = data;
            else
                out.freq(f).pt(p).sp_or_not(sp).data = data(:,:,:,f);
            end
        end
        
    end
end

end