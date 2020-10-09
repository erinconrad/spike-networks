function count_good_spikes

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
bad_spikes_folder = [results_folder,'bad_spikes/'];

listing = dir([bad_spikes_folder,'*.mat']);

all_names = {};
sp_counts = [];
not_counts = [];
for i = 1:length(listing)
    filename = listing(i).name;
    name_sp = split(filename,'_');
    name = name_sp{1};
    
    
    
    [a,b] = ismember(name,all_names);
    if a == 1
        pt_idx = b;
    else
        all_names = [all_names;name];
        pt_idx = length(all_names);
    end
    
    % Load the file
    bad = load([bad_spikes_folder,filename]);
    bad = bad.bad;
    
    count = 0;
    for s = 1:length(bad.spike)
        if bad.spike(s).bad_spike == 0
            count = count + 1;
        end
    end
    
    if contains(filename,'not') == 1
        not_counts = [not_counts;count];
    else
        sp_counts = [sp_counts;count];
    end
end


table(all_names,sp_counts,not_counts,'variablenames',...
    {'name','spike','notspike'})

end