function ers_change(time_window)

%% Parameters

freq_bands = [0 4;... %delta
    4 8;...%theta
    8 12;...% alpha
    12 24;... %beta
    30 40;... % low gamma
    96 106;... % high gamma
    106 256;... %ultra-high
    0 256;... %broadband    
    ]; 
freq_names = {'delta','theta','alpha','beta','low_gamma',...
    'high_gamma','ultra_high','broadband'};
n_f = size(freq_bands,1);
time_text = sprintf('%1.1f/',time_window);

%% Get file locations, load spike times and pt structure
locations = spike_network_files;
main_folder = locations.main_folder;
eeg_folder = [main_folder,'results/eeg_data/'];
results_folder = [main_folder,'results/'];
data_folder = [main_folder,'data/'];
script_folder = locations.script_folder;
addpath(genpath(script_folder));
pt_file = [data_folder,'spike_structures/pt.mat'];
ers_folder = [results_folder,'ers/',time_text,'/'];


end