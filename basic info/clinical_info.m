%% Get file locations, load spike times and pt structure
locations = spike_network_files;
main_folder = locations.main_folder;
results_folder = [main_folder,'results/'];
script_folder = locations.script_folder;
addpath(genpath(script_folder));
pt_folder = [main_folder,'data/spike_structures/'];
power_folder = [results_folder,'power/manual/0.1/'];

% load pt file
pt = load([pt_folder,'pt.mat']);
pt = pt.pt;


% load power file (because it has all the included patients)
power = load([power_folder,'sig_dev.mat']);
power = power.sig_dev;
npts = length(power);

% Get pt indices from pt file
pt_indices = [];
for i = 1:length(power)
    name = power(i).name;
    for p = 1:length(pt)
        if strcmp(name,pt(p).name)
            pt_indices = [pt_indices;p];
        end
    end
end

if length(pt_indices) ~= npts, error('what'); end

% Get some clinical info
sex = {};
age = [];

for i = 1:length(pt_indices)
    p = pt_indices(i);
    sex = [sex;pt(p).clinical.sex];
    age = [age;str2num(pt(p).clinical.ageSurgery)];
end

fprintf('\nThere were %d females.\n',sum(strcmp(sex,'F')));
fprintf('\nThe mean age was %1.1f (range %d-%d)\n',mean(age),min(age),max(age));
