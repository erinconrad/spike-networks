function out = identify_bad_spikes


%% Get file locations, load spike times and pt structure
locations = spike_network_files;
main_folder = locations.main_folder;
results_folder = [main_folder,'results/'];
bad_folder = [results_folder,'bad_spikes/'];


listing = dir([bad_folder,'*.mat']);
%% Get patient names
all_names = {};
for i = 1:length(listing)
    fname = listing(i).name;
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
    pt_listing = dir([bad_folder,'*',name,'*']);
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
        bad = load([bad_folder,sp_or_not_listing{sp}]);
        bad = bad.bad;
        
        %% Get bad spikes
        bad_spikes = zeros(length(bad.spike),1);
        for s = 1:length(bad.spike)
            if bad.spike(s).bad_spike == 1
                bad_spikes(s) = 1;
            end
        end
        bad_spikes = logical(bad_spikes);
        
        %% Restructure
        out.pt(p).name = name;
        out.pt(p).sp_or_not(sp).bad = bad_spikes;
        
    end
    
end



end