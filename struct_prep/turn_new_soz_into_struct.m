function out_labels = turn_new_soz_into_struct(name)

%% Get file locations, load spike times and pt structure
locations = spike_network_files;
main_folder = locations.main_folder;
data_folder = [main_folder,'data/'];
results_folder = [main_folder,'results/'];
data_folder = [main_folder,'data/'];
eeg_folder = [main_folder,'results/eeg_data/'];
script_folder = locations.script_folder;
addpath(genpath(script_folder));

soz_folder = [data_folder,'new_soz/'];

filename = [soz_folder,'SOZ.xlsx'];

T = readtable(filename);
soz_labels = T.(name);

out_labels = {};
% remove empty
for i = 1:length(soz_labels)
    if ~isempty(soz_labels{i})
        curr_label = soz_labels{i};
        
        %% Remove leading zero
        % Get number
        numIdx = regexp(curr_label,'\d');
        
        % get first
        numIdx = numIdx(1);
        
        % see if it's zero
        if strcmp(curr_label(numIdx),'0')
            
            % remove it
            curr_label(numIdx) = [];
            
        end
        
        out_labels = [out_labels;curr_label];
    end
end

end