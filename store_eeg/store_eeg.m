function store_eeg(times,pt)

%% Parameters
surround_time = 5;


%% Get locations
locations = spike_network_files;
script_folder = locations.script_folder;
main_folder = locations.main_folder;
pwname = locations.pwfile;
addpath(genpath(script_folder));

results_folder =  [main_folder,'results/'];
if exist(results_folder,'dir') == 0
    mkdir(results_folder);
end


for whichPt = 1:length(times)
    
    name = times(whichPt).name;
    
    if isempty(name) == 1
        continue;
    else
        fprintf('Doing %s...\n',name);
    end
    
    % Make patient folder
    pt_folder = [results_folder,name,'/'];
    if exist(pt_folder,'dir') == 0
        mkdir(pt_folder);
    end
    
    % ensure it is the same patient as in pt
    if strcmp(name,pt(whichPt).name) == 0
        error('what\n');
    end
    
    
    
    % get ieeg name and sampling rate
    ieeg_name = pt(whichPt).ieeg_name;
    fs = pt(whichPt).fs;
    
    % get which channels to do: this tells me which ieeg channels to take
    % and in what order so that they should perfectly align with the
    % electrode designations in pt.new_elecs
    chs = pt(whichPt).new_elecs.ch_order;
    
    %% Do a test download to make sure the channels line up perfectly
    data = download_eeg(ieeg_name,[],pwname,1);
    if data.fs ~= fs, error('what\n'); end
    % Save basic info about the patient
    save([pt_folder,'basic_info.mat'],'data');
    
    %% Loop through spike times
    for do_spike = [1] % hold off on not-a-spike times
        
        if do_spike == 1
            s_times = times(whichPt).spike_times;
            fname = 'spikes.mat';
        elseif do_spike == 2
            s_times = times(whichPt).not_spike_times;
            fname = 'not_spikes.mat';
        end
        
        for t = 1:length(s_times)
                         
            % which times
            which_times = [s_times(t)-surround_time,s_times(t)+surround_time];
            
            % get indices
            indices = round(which_times(1)*fs):round(which_times(2)*fs);

            % Get the data
            data = download_eeg(ieeg_name,indices,pwname,0); 
            
            spike(t).values = data.values;
            spike(t).time = s_times(t);
            spike(t).label = times(whichPt).spike_labels{t};
            
        end
        
    end

    % Save the data in a .mat file
    save([pt_folder,fname],'spike');
    
    
    
end


end