function spike_stats(whichPts)

%% Parameters
% 1 = alpha/theta; 2 = beta, 3 = low gamma, 4 = high gamma, 5 = ultra high, 6 = broadband
freq_text = {'alpha/theta','beta','low\ngamma','high\ngamma','ultra high\ngamma','broadband'};
%freq_text = {'alpha/theta'};
n_f = length(freq_text);
n_times = 11;

%% Get file locations, load spike times and pt structure
locations = spike_network_files;
main_folder = locations.main_folder;
results_folder = [main_folder,'results/'];
data_folder = [main_folder,'data/'];
script_folder = locations.script_folder;
addpath(genpath(script_folder));
spike_times_file = [data_folder,'spike_times/times.mat'];
pt_file = [data_folder,'spike_structures/pt.mat'];

times = load(spike_times_file); % will result in a structure called "out"
times = times.out;
pt = load(pt_file); % will create a structure called "pt"
pt = pt.pt;

if isempty(whichPts) == 1
    whichPts = 1:length(times);
end

for whichPt = whichPts
    
    %% Prep patient
    
    % Skip it if name is empty
    if isempty(times(whichPt).name) == 1, continue; end
    name = times(whichPt).name;
    fprintf('\nDoing %s\n',name);
    
    pt_folder = [results_folder,name,'/'];
    stats_folder = [pt_folder,'stats/'];
    
    % Load the stats file
    out = load([stats_folder,'stats.mat']);
    stats = out.out;
    
    % Get variables
    bin_z = stats.signal.bin_z;
    ec = stats.network.ec;
    
    % Remove spike second
    bin_z_nospike = bin_z(:,[1:5,7:11]);
    ec_nospike = ec(:,:,[1:5,7:11]);
    
    % Just look at alpha/theta
    ec_nospike_at = squeeze(ec(1,:,:));
    
    %% Do stats
    
    % Do a repeated measures ANOVA
    t = table(bin_z_nospike);
    rm = fitrm(t);
    
end

end