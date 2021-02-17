function [metrics,is_spike_soz,pre_spike,pt_rise,is_spike_depth] = remove_bad_spikes(metrics,is_spike_soz,pre_spike,pt_rise,is_spike_depth)

old_metrics = metrics;

%% Get file locations, load spike times and pt structure
locations = spike_network_files;
main_folder = locations.main_folder;
results_folder = [main_folder,'results/'];
data_folder = [main_folder,'data/'];
eeg_folder = [main_folder,'results/eeg_data/'];
script_folder = locations.script_folder;
addpath(genpath(script_folder));

bct_folder = locations.BCT;
addpath(genpath(bct_folder));
out_folder = [results_folder,'plots/'];
sig_dev_folder = [results_folder,'signal_deviation/manual/'];
perm_folder = [results_folder,'perm_stats/'];
ers_folder = [results_folder,'ers/'];
ns_folder = [results_folder,'metrics/manual/'];
sp_diff_folder = [results_folder,'net_diff_stats/'];
bad_folder = [results_folder,'bad_spikes/'];

% Loop through patients
all_names = {};
listing = dir([bad_folder,'*.mat']);
for i = 1:length(listing)
    filename = listing(i).name;
    pt_name = strsplit(filename,'_');
    pt_name = pt_name{1};
    
    [a,b] = ismember(pt_name,all_names);
    if a == 1
        pt_idx = b;
    else
        all_names = [all_names;pt_name];
        pt_idx = length(all_names);
    end
    
    % Load pt file
    bad = load([bad_folder,filename]);
    bad = bad.bad;
    
    % Get bad spikes
    bad_spikes = zeros(length(bad.spike),1);
    for s = 1:length(bad.spike)
        if bad.spike(s).bad_spike == 1
            bad_spikes(s) = 1;
        end
    end
    bad_spikes = logical(bad_spikes);
    
    if ~contains(filename,'not')
        pre_spike(pt_idx).windows.manual_before_rise(bad_spikes,:) = [];
    end
    
    % Loop through metrics
    for n = 1:length(metrics)
        for t = 1:length(metrics(n).time)
            for f = 1:length(metrics(n).time(t).freq)
                fn = fieldnames(metrics(n).time(t).freq(f));
                for ff = 1:length(fn)
                    met = fn{ff};
                    if strcmp(met,'name'), continue; end
                    
                    
                    
                    % Remove the bad spikes
                    if contains(filename,'not')
                        if size(bad_spikes,1) > size(metrics(n).time(t).freq(f).(met).pt(pt_idx).not.data,1)
                            bad_spikes = bad_spikes(1:size(metrics(n).time(t).freq(f).(met).pt(pt_idx).not.data,1),:);
                        end
                        metrics(n).time(t).freq(f).(met).pt(pt_idx).not.data(bad_spikes,:) = [];
                    else
                        if size(bad_spikes,1) > size(metrics(n).time(t).freq(f).(met).pt(pt_idx).spike.data,1)
                            bad_spikes = bad_spikes(1:size(metrics(n).time(t).freq(f).(met).pt(pt_idx).spike.data,1),:);
                        end
                        metrics(n).time(t).freq(f).(met).pt(pt_idx).spike.data(bad_spikes,:) = [];
                        metrics(n).time(t).freq(f).(met).pt(pt_idx).first(bad_spikes,:) = [];
                        metrics(n).time(t).freq(f).(met).pt(pt_idx).other(bad_spikes,:) = [];
                    end
                end
            end
        end
    end
    
    % Also remove from soz struct
    if ~contains(filename,'not')
        is_spike_soz(pt_idx).is_soz(bad_spikes) = [];
    end
    
    % Pt rise
    if ~contains(filename,'not')
        pt_rise(pt_idx).both(bad_spikes,:) = [];
    end
    
    if ~contains(filename,'not')
        is_spike_depth(pt_idx).is_depth(bad_spikes) = [];
        
    
end



end