function show_involved_chs(p)


%% Parameters
surround_time = 2;

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
power_folder = [results_folder,'power/manual/0.1/'];

%% Load power file (because this gets names from index)
power = load([power_folder,'sig_dev.mat']);
power = power.sig_dev;
% Get name
name = power(p).name;

% load eeg data
filename = [name,'_eeg.mat'];
spike = load([eeg_folder,filename]);
spike = spike.spike;
values = spike(1).data;
nchs = size(values,2);
fs = spike(1).fs;

% loop through spikes
for s = 1:length(spike)

    % get eeg data
    data = spike(s).data; % ntimes x nch

    % pre process
    data = pre_processing(data,do_car,pre_whiten,do_notch,fs);

    % get involved chs
    is_sp_ch = spike(s).involved;
    
    % get order of chs
    ordered_chs = spike(s).ordered_chs(:,1);
    
    biggest_dev_manual = spike(s).biggest_dev;
    
    
    figure
    set(gcf,'position',[100 100 1000 500])
    offset = 0;
    ch_offsets = zeros(sum(is_sp_ch),1);
    ch_bl = zeros(sum(is_sp_ch),1);
    for ich = 1:size(spike(s).ordered_chs,1)
        if ordered_chs(ich) == biggest_dev_manual
            plot(linspace(-3,3,size(data,1)),data(:,ordered_chs(ich))+offset,'linewidth',3)
        else
            plot(linspace(-3,3,size(data,1)),data(:,ordered_chs(ich))+offset)
        end
        ch_offsets(ich) = offset;
        ch_bl(ich) = offset + median(data(:,ordered_chs(ich)));
        hold on
        text(surround_time+0.05,ch_bl(ich),sprintf('%s',spike(s).chLabels{ordered_chs(ich)}))        
        if ich<sum(is_sp_ch)
            offset = offset + max(data(:,ordered_chs(ich))) - min(data(:,ordered_chs(ich+1)));
        end
            
    end
    xlim([-surround_time,surround_time]);
    pause
    close(gcf)
    

end

