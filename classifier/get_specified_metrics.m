function [stats,is_spike_soz,is_spike_depth] = get_specified_metrics(windows,pre_spike,met)

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
ns_folder = [results_folder,'metrics/manual/coherence/'];
sp_diff_folder = [results_folder,'net_diff_stats/coherence/'];
biggest_dev_folder = [results_folder,'biggest_dev/'];
seq_folder = [results_folder,'seq_data/'];


%% load pt file
pt = load([data_folder,'spike_structures/pt.mat']);
pt = pt.pt;

%% Specify metric folders
time_name = sprintf('%1.1f/',windows);

if contains(met,'ers')
    main_folder = ers_folder;
elseif strcmp(met,'F')
    main_folder =  sp_diff_folder;
elseif strcmp(met,'ns_inv') || strcmp(met,'ge') || strcmp(met,'ns_big') || ...
                strcmp(met,'ns_avg') || strcmp(met,'ns_auto')
    main_folder = ns_folder;
elseif contains(met,'sd')
    main_folder = ers_folder;
end

time_folder = [main_folder,time_name];

pt_listing = dir([time_folder,'*.mat']);

% for absolute power, fill up the info from the pre_spike structure that we
% already got
if contains(met,'sd')
    
        
    for f = 1:3


        %times = round((ers.index_windows(:,1)/ers.fs-3)*1e2)/(1e2);
        for p = 1:length(pre_spike)
            stats.time.freq(f).(met).pt(p).times = pre_spike(p).windows.cons_windows;
            %{
            if strcmp(wpr,'cons')
                stats(n).time(t).freq(f).sd.pt(p).times = pre_spike(p).windows(t).cons_windows;
            else
                stats(n).time(t).freq(f).sd.pt(p).times = pre_spike(p).windows(t).all_windows;
            end
            %}
            stats.time.freq(f).(met).pt(p).name = pre_spike(p).name;
            stats.time.freq(f).(met).pt(p).spike.data = pre_spike(p).windows.dev.spike;
            stats.time.freq(f).(met).pt(p).not.data = pre_spike(p).windows.dev.not;
            stats.time.freq(f).(met).pt(p).first = pre_spike(p).windows.dev.first;
            stats.time.freq(f).(met).pt(p).other = pre_spike(p).windows.dev.other;
        end
    end
    
end

all_names = {};
        % loop through pts
for i = 1:length(pt_listing)

    fname = pt_listing(i).name;
    pt_name = strsplit(fname,'_');
    pt_name = pt_name{1};


    [a,b] = ismember(pt_name,all_names);
    if a == 1
        pt_idx = b;
    else
        all_names = [all_names;pt_name];
        pt_idx = length(all_names);
    end

    % load pt file
    load([time_folder,fname]); % this will produce a structure of different names depending on metric
    
    % Load spike file
    %{
    spike = load([ers_folder,time_name,sprintf('%s_ers.mat',pt_name)]);
    spike = spike.ers; % ers file
    %}
    spike = load([eeg_folder,pt_name,'_eeg.mat']);
    
    if contains(fname,'not') == 0
        manual_big = load([biggest_dev_folder,pt_name,'_rise.mat']);
        manual_big = manual_big.early;
        if length(manual_big.spike) ~= length(spike.spike)
            error('what');
        end
        
        seq = load([seq_folder,pt_name,'_seq.mat']);
        seq = seq.seq;
        if length(seq) ~= length(spike.spike)
            error('what');
        end
    end

    % Get number of frequencies
    if contains(met,'ers')
        nfreq = size(ers.freq_bands,1);
        for f = 1:nfreq
            stats(1).time(1).freq(f).name = ers.freq_names{f};

        end
    elseif strcmp(met,'ns_inv') || strcmp(met,'ns_auto') || strcmp(met,'ns_big') || ...
                strcmp(met,'ns_avg') || strcmp(met,'ge')
        nfreq = length(metrics.freq);
        for f = 1:nfreq
            stats(1).time(1).freq(f).name = metrics.freq(f).name;
        end
    elseif contains(met,'sd')
        nfreq = size(ers.freq_bands,1);
        for f = 1:nfreq
            stats(1).time(1).freq(f).name = ers.freq_names{f};

        end
    else

        nfreq = length(sim);
        for f = 1:nfreq
            stats(1).time(1).freq(f).name = sim.freq(f).name;
        end
    end
    stats(1).time(1).time_window = windows;
    
    %% Get soz channels
    if contains(fname,'not') == 0
        soz_chs = get_soz_chs(pt,pt_name);
    elseif contains(fname,'not') == 1
    else
        is_spike_soz = [];
    end

    for f = 1:nfreq
        if strcmp(met,'ns_big') || strcmp(met,'ns_avg') || strcmp(met,'ns_auto') || strcmp(met,'ge')

            sim = metrics;
            ns = sim.freq(f).ns.data;
            ns_all = sim.freq(f).ns_all.data;
            ns_auto = sim.freq(f).ns_auto.data;
            ge = sim.freq(f).ge.data;
            

            % Avg ns_all across channels
            ns_avg = squeeze(mean(ns_all,3));

            times = round((sim.index_windows(:,1)/sim.fs-3)*1e2)/(1e2);
            
            stats(1).time(1).freq(f).ns_avg.pt(pt_idx).name = pt_name;
            stats(1).time(1).freq(f).ns_big.pt(pt_idx).name = pt_name;

            % spike vs not a spike
            if contains(fname,'not') == 1
                stats(1).time(1).freq(f).ns_avg.pt(pt_idx).not.data(:,:) = ns_avg;
                stats(1).time(1).freq(f).ns_big.pt(pt_idx).not.data(:,:) = ns;
                stats(1).time(1).freq(f).ns_auto.pt(pt_idx).not.data(:,:) = ns_auto;
                stats(1).time(1).freq(f).ge.pt(pt_idx).not.data(:,:) = ge;

            else
                ns_first = sim.freq(f).ns_first.data;
                ns_other = sim.freq(f).ns_other.data;
                
                stats(1).time(1).freq(f).ns_avg.pt(pt_idx).spike.data(:,:) = ns_avg;
                stats(1).time(1).freq(f).ns_big.pt(pt_idx).spike.data(:,:) = ns;
                stats(1).time(1).freq(f).ns_auto.pt(pt_idx).spike.data(:,:) = ns_auto;
                stats(1).time(1).freq(f).ge.pt(pt_idx).spike.data(:,:) = ge;
                
                stats(1).time(1).freq(f).ns_avg.pt(pt_idx).first = ns_first;
                stats(1).time(1).freq(f).ns_avg.pt(pt_idx).other = ns_other;
                
                stats(1).time(1).freq(f).ns_big.pt(pt_idx).first = ns_first;
                stats(1).time(1).freq(f).ns_big.pt(pt_idx).other = ns_other;
                
                stats(1).time(1).freq(f).ns_auto.pt(pt_idx).first = ns_first;
                stats(1).time(1).freq(f).ns_auto.pt(pt_idx).other = ns_other;
                
                stats(1).time(1).freq(f).ge.pt(pt_idx).first = ns_first;
                stats(1).time(1).freq(f).ge.pt(pt_idx).other = ns_other;
            end

            stats(1).time(1).freq(f).ns_big.name = 'Node strength (spike channel)';
            stats(1).time(1).freq(f).ns_avg.name = 'Node strength (average)';
          

            stats(1).time(1).freq(f).ns_big.pt(pt_idx).times = times;
            stats(1).time(1).freq(f).ns_avg.pt(pt_idx).times = times;
            stats(1).time(1).freq(f).ns_auto.pt(pt_idx).times = times;
            stats(1).time(1).freq(f).ge.pt(pt_idx).times = times;
        end

        % Add ERS stuff
        if contains(met,'ers')

            stats(1).time(1).freq(f).(met).pt(pt_idx).index_windows = ers.index_windows;
            times = round((ers.index_windows(:,1)/ers.fs-3)*1e2)/(1e2);
            stats(1).time(1).freq(f).(met).pt(pt_idx).times = times;
            stats(1).time(1).freq(f).(met).pt(pt_idx).name = ers.name;
            
            all_ers(f).data = zeros(length(ers.spike),size(ers.spike(1).ers,1)); % n spikes x n times
            all_ers(f).first = nan(length(ers.spike),size(ers.spike(1).ers,1));
            all_ers(f).other = nan(length(ers.spike),size(ers.spike(1).ers,1));
            for s = 1:length(ers.spike)
                if ~isempty(ers.spike(s).ers_first)
                    all_ers(f).first(s,:) = ers.spike(s).ers_first(:,f);
                    all_ers(f).other(s,:) = ers.spike(s).ers_others(:,f);
                end
            end
            if strcmp(met,'ers')
                
                for s = 1:length(ers.spike)
                    all_ers(f).data(s,:) = ers.spike(s).ers(:,f);
                    
                end
            elseif strcmp(met,'ers_auto')
                for s = 1:length(ers.spike)
                    all_ers(f).data(s,:) = ers.spike(s).ers_auto(:,f);
                end
           
            end
            
           

            if contains(fname,'not') == 1
                stats(1).time(1).freq(f).(met).pt(pt_idx).not.data = all_ers(f).data;
            else
                stats(1).time(1).freq(f).(met).pt(pt_idx).spike.data = all_ers(f).data;
                stats(1).time(1).freq(f).(met).pt(pt_idx).first = all_ers(f).first;
                stats(1).time(1).freq(f).(met).pt(pt_idx).other = all_ers(f).other;
            end
        end

        

    end
    
    %% Add soz info
    if contains(fname,'not') == 0
        %is_soz_pt = zeros(length(ers.spike),1);
        is_soz_pt = zeros(length(spike.spike),1);
        is_depth = zeros(length(spike.spike),1);
        for s = 1:length(spike.spike)
            %
            if contains(met,'auto')
                biggest_dev = spike.spike(s).biggest_dev;
            else
                biggest_dev = manual_big.spike(s).dev_ch;
            end


            
            if ismember(biggest_dev,soz_chs)
                is_soz_pt(s) = 1;
            else
                is_soz_pt(s) = 0;
            end

            is_depth(s) = is_it_depth(pt,biggest_dev,pt_name);
        end
        is_spike_soz(pt_idx).is_soz = is_soz_pt;
        is_spike_depth(pt_idx).is_depth = is_depth;
    end


end

    





end