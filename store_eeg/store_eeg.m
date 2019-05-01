function store_eeg(times,pt)

%% Parameters
surround_time = 7;


%% Get locations
locations = spike_network_files;
script_folder = locations.script_folder;
main_folder = locations.main_folder;
pwname = locations.pwfile;
addpath(genpath(script_folder));

if isempty(locations.ieeg_folder) == 0
    addpath(genpath(locations.ieeg_folder));
end

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
    
    
    % How many spikes does this patient have
    n_spikes = sum(~isnan(times(12).spike_times));
    n_files = floor(n_spikes/100)+1;
    
    % See what files are already present, skip if all done
    listing = dir([pt_folder,'*.mat']);
    if length(listing) == n_files + 1 % 1 extra for the info file
        fprintf('Already did %s, skipping...\n',name);
        continue
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
    if isempty(listing) == 1
        data = download_eeg(ieeg_name,[],pwname,1);
        if data.fs ~= fs, error('what\n'); end
        % Save basic info about the patient
        save([pt_folder,'basic_info.mat'],'data');
    end
    
    %% Loop through spike times
    for do_spike = [1] % hold off on not-a-spike times
        
        if do_spike == 1
            s_times = times(whichPt).spike_times;
            fname = 'spikes';
        elseif do_spike == 2
            s_times = times(whichPt).not_spike_times;
            fname = 'not_spikes';
        end
        
        if length(listing) > 1
            n_spikes_done = 100*(length(listing)-1);
            f_index = length(listing)-1;
        else
            n_spikes_done = 0;
            f_index = 0;
        end
        
        % start with first we haven't done
        count = 0;
        for t = n_spikes_done+1:length(s_times) 
            
            count = count + 1;
            
            if isnan(s_times(t)) == 0
           
                         
                % which times
                which_times = [s_times(t)-surround_time,s_times(t)+surround_time];

                % get indices
                indices = round(which_times(1)*fs):round(which_times(2)*fs);

                % Get the data
                data = download_eeg(ieeg_name,indices,pwname,0); 

                spike(count).values = data.values;
                spike(count).time = s_times(t);
                spike(count).label = times(whichPt).spike_labels{t};
                spike(count).which = t;
                temp_cell = times(whichPt).seq_labels(t,:);
                spike(count).seq_labels = temp_cell(~cellfun(@isempty,temp_cell));
                
            end
            
            % save when we get to 100 and restart the count
            if count == 100
                fprintf('Did %d spikes of %s...\n',t,name);
                count = 0;
                f_index = f_index + 1;
                save([pt_folder,fname,'_',sprintf('%d',f_index),'.mat'],'spike');
            end
            
        end
        
    end

    
    
end


end