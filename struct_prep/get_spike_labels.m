function get_spike_labels

%% Get locations
locations = spike_network_files;
script_folder = locations.script_folder;
main_folder = locations.main_folder;
pwname = locations.pwfile;
data_folder = [main_folder,'data/'];
pt_folder = [data_folder,'spike_structures/'];
results_folder = [main_folder,'results/'];
addpath(genpath(script_folder));
eeg_folder = [results_folder,'eeg_data/'];
manual_dev_folder = [results_folder,'biggest_dev/'];
power_folder = [results_folder,'power/manual/0.1/'];
out_folder = [results_folder,'ied_labels/'];
if exist(out_folder,'dir') == 0
    mkdir(out_folder);
end

%% Load power file (because this gets names from index)
power = load([power_folder,'sig_dev.mat']);
power = power.sig_dev;

all_pt_chs = cell(length(power),1);
all_pt_names = cell(length(power),1);
max_num_chs = 0;
for p = 1:length(power)
    % Get name
    name = power(p).name;

    %% Load eeg file
    eeg_file = [eeg_folder,name,'_eeg.mat'];
    auto = load(eeg_file);
    auto = auto.spike;

    % Get channel labels for auto biggest dev
    auto_chs = cell(length(auto),1);
    for s = 1:length(auto)
        auto_chs{s} = auto(s).chLabels{auto(s).biggest_dev};
    end
    auto_chs = unique(auto_chs);

    
    %% Load manual biggest dev
    manual_file = [manual_dev_folder,name,'_rise.mat'];
    manual = load(manual_file);
    manual = manual.early;

    % Get channel labels for manual biggest dev
    manual_chs = cell(length(manual.spike),1);
    for s = 1:length(manual.spike)
        manual_chs{s} = manual.spike(s).dev_ch_label;
    end

    manual_chs = unique(manual_chs);

    if 0
    %% Display channel names
    fprintf('\nAutomatically determined IED peak electrodes:\n');
    for i = 1:length(auto_chs)
        fprintf('\n%s',auto_chs{i});
    end

    fprintf('\nManually determined IED peak electrodes:\n');
    for i = 1:length(manual_chs)
        fprintf('\n%s',manual_chs{i});
    end
    fprintf('\n');
    end
    
    all_chs = unique([manual_chs;auto_chs]);
    all_pt_chs{p} = all_chs;
    all_pt_names{p} = name;
    
    if length(all_chs) > max_num_chs
        max_num_chs = length(all_chs);
    end
    
end

% Pad channel labels
for p = 1:length(all_pt_chs)
    ch_num = length(all_pt_chs{p});
    too_short = max_num_chs - ch_num;
    pad = cell(too_short,1);
    for i = 1:length(pad)
        pad{i} = '';
    end
    all_pt_chs{p} = [all_pt_chs{p};pad];
end

new_cell = cell(max_num_chs,length(all_pt_chs));
for i = 1:size(new_cell,2)
    new_cell(:,i) = (all_pt_chs{i});
end

new_new_cell = cell(max_num_chs,length(all_pt_chs)*2);
for i = 1:size(new_new_cell,2)
    if mod(i,2) == 1
        new_new_cell(:,i) = (all_pt_chs{floor(i/2)+1});
    else
        new_new_cell(:,i) = num2cell(zeros(max_num_chs,1));
    end
end

out_names = cell(1,length(all_pt_chs)*2);
for i = 1:size(out_names,2)
    if mod(i,2) == 1
        out_names{i} = all_pt_names{floor(i/2)+1}; 
    else
        out_names{i} = [all_pt_names{floor(i/2)},'ana']; 
    end
end

T = cell2table(new_new_cell,'VariableNames',out_names);
writetable(T,[out_folder,'labels.csv'])


end