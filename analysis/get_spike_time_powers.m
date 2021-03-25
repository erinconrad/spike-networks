function out = get_spike_time_powers(metrics)

%% Get locations
locations = spike_network_files;
main_folder = locations.main_folder;
results_folder = [main_folder,'results/spike_time_power/'];

for p = 1:length(metrics.freq(1).pt)
    name = metrics.freq(1).pt(p).name;
    out(p).name = name;
    
    % load
    power = load([results_folder,name,'_spike_time_power.mat']);
    
    out(p).spike_powers = power.spike_powers;
end

end