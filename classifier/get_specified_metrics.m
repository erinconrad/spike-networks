function [stats,is_spike_soz] = get_specified_metrics(windows,pre_spike,met)

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


%% Specify metric folders
time_name = sprintf('%1.1f/',windows);

if contains(met,'ers')
    main_folder = ers_folder;
elseif strcmp(met,'F')
    main_folder =  sp_diff_folder;
elseif strcmp(met,'ns_inv') || strcmp(met,'ge') || strcmp(met,'ns_big') || ...
                strcmp(met,'ns_avg') || strcmp(met,'trans')
    main_folder = ns_folder;
elseif contains(met,'sd')
    main_folder = ers_folder;
end

time_folder = [main_folder,time_name];

pt_listing = dir([time_folder,'*.mat']);


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
    spike = load([ers_folder,time_name,sprintf('%s_ers.mat',pt_name)]);
    spike = spike.ers; % ers file
    
    if contains(fname,'not') == 0
        manual_big = load([biggest_dev_folder,pt_name,'_rise.mat']);
        manual_big = manual_big.early;
        if length(manual_big.spike) ~= length(spike.spike)
            error('what');
        end
    end

    % Get number of frequencies
    if contains(met,'ers')
        nfreq = size(ers.freq_bands,1);
        for f = 1:nfreq
            stats(1).time(1).freq(f).name = ers.freq_names{f};

        end
    elseif strcmp(met,'F')
        nfreq = length(sim);
        for f = 1:nfreq
            stats(1).time(1).freq(f).name = sim.name;
        end
    elseif strcmp(met,'ns_inv') || strcmp(met,'ns_auto') || strcmp(met,'ns_big') || ...
                strcmp(met,'ns_avg')
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
        soz_chs = get_soz_chs(pt_name);
    elseif contains(fname,'not') == 1
    else
        is_spike_soz = [];
    end

    for f = 1:nfreq
        if strcmp(met,'ns_big') || strcmp(met,'ns_avg') || strcmp(met,'ns_auto')

            sim = metrics;
            ns = sim.freq(f).ns.data;
            ns_all = sim.freq(f).ns_all.data;
            ns_auto = sim.freq(f).ns_auto.data;

            % Avg ns_all across channels
            ns_avg = squeeze(mean(ns_all,3));

            times = round((sim.index_windows(:,1)/sim.fs-3)*1e2)/(1e2);
            
            stats(1).time(1).freq(f).ns_avg.pt(pt_idx).name = pt_name;
            stats(1).time(1).freq(f).ns_big.pt(pt_idx).name = pt_name;

            % spike vs not a spike
            if contains(fname,'not') == 1
                stats(1).time(1).freq(f).ns_avg.pt(pt_idx).not.data(:,:) = ns_avg;
                stats(1).time(1).freq(f).ns_big.pt(pt_idx).not.data(:,:) = ns;
                stats(1).time(1).freq(f).ns_auto.pt(pt_idx).not.data(:,:) = ns;

            else
                stats(1).time(1).freq(f).ns_avg.pt(pt_idx).spike.data(:,:) = ns_avg;
                stats(1).time(1).freq(f).ns_big.pt(pt_idx).spike.data(:,:) = ns;
                stats(1).time(1).freq(f).ns_auto.pt(pt_idx).spike.data(:,:) = ns;
            end

            stats(1).time(1).freq(f).ns_big.name = 'Node strength (spike channel)';
            stats(1).time(1).freq(f).ns_avg.name = 'Node strength (average)';
          

            stats(1).time(1).freq(f).ns_big.pt(pt_idx).times = times;
            stats(1).time(1).freq(f).ns_avg.pt(pt_idx).times = times;
            
        end

        % Add ERS stuff
        if contains(met,'ers')

            stats(1).time(1).freq(f).(met).pt(pt_idx).index_windows = ers.index_windows;
            times = round((ers.index_windows(:,1)/ers.fs-3)*1e2)/(1e2);
            stats(1).time(1).freq(f).(met).pt(pt_idx).times = times;
            stats(1).time(1).freq(f).(met).pt(pt_idx).name = ers.name;
            
            all_ers(f).data = zeros(length(ers.spike),size(ers.spike(1).ers,1)); % n spikes x n times
            if strcmp(met,'ers')
                
                for s = 1:length(ers.spike)
                    all_ers(f).data(s,:) = ers.spike(s).ers(:,f);
                end
            elseif strcmp(met,'ers_auto')
                for s = 1:length(ers.spike)
                    all_ers(f).data(s,:) = ers.spike(s).ers_auto(:,f);
                end
            elseif strcmp(met,'rel_big_dev_ers')
                for s = 1:length(ers.spike)
                    all_ers(f).data(s,:) = ers.spike(s).ers_all(:,f,ers.spike(s).biggest_dev)./...
                        nanmean(ers.spike(s).ers_all(:,f,:),3);
                end
            elseif strcmp(met,'ers_avg')
                for s = 1:length(ers.spike)
                    all_ers(f).data(s,:) = nanmean(ers.spike(s).ers_all(:,f,:),3);
                        
                end
            elseif strcmp(met,'ers_soz')
                for s = 1:length(ers.spike)
                    all_ers(f).data(s,:) = nanmean(ers.spike(s).ers_all(:,f,soz_chs),3);    
                end
            elseif strcmp(met,'ers_only_keep_soz')
                for s = 1:length(ers.spike)
                    biggest_dev = ers.spike(s).biggest_dev;
                    if ismember(biggest_dev,soz_chs)
                        all_ers(f).data(s,:) = nanmean(ers.spike(s).ers_all(:,f,soz_chs),3);    
                    else
                        all_ers(f).data(s,:) = nan;  
                    end
                end
            end
            
           

            if contains(fname,'not') == 1
                stats(1).time(1).freq(f).(met).pt(pt_idx).not.data = all_ers(f).data;
            else
                stats(1).time(1).freq(f).(met).pt(pt_idx).spike.data = all_ers(f).data;
            end
        end

        % Add high gamma power stuff
        if strcmp(met,'hg')
            stats(1).time(1).freq(f).hg.pt(pt_idx).index_windows = hg.index_windows;
            times = round((hg.index_windows(:,1)/hg.fs-3)*1e2)/(1e2);
            stats(1).time(1).freq(f).hg.pt(pt_idx).times = times;
            stats(1).time(1).freq(f).hg.pt(pt_idx).name = ers.name;

            if contains(fname,'not') == 1
                stats(1).time(1).freq(f).hg.pt(pt_idx).not.data = squeeze(mean(hg.powers,4));
            else
                stats(1).time(1).freq(f).hg.pt(pt_idx).spike.data = squeeze(mean(hg.powers,4));
            end
        end

        % Add spike diff stuff
        if strcmp(met,'F')
            spd.sim = sim;
            if f>length(spd.sim)
                fprintf('\nWarning, non-aligning frequencies\n');
                continue;
            end
            stats(1).time(1).freq(f).F.pt(pt_idx).index_windows = spd.sim(f).index_windows;
            times = round((spd.sim(f).index_windows(:,1)/spd.sim(f).fs-3)*1e2)/(1e2);
            stats(1).time(1).freq(f).score.pt(pt_idx).index_windows = spd.sim(f).index_windows;
            stats(1).time(1).freq(f).score.pt(pt_idx).times = times;
            stats(1).time(1).freq(f).F.pt(pt_idx).times = times;

            stats(1).time(1).freq(f).F.pt(pt_idx).name = spd.sim(f).pt_name;
            stats(1).time(1).freq(f).score.pt(pt_idx).name = spd.sim(f).pt_name;
            % convert score to a 2 dimensional vector
            score = spd.sim(f).score;
            time_idx = spd.sim(f).time_idx;
            num_time_idx = length(unique(time_idx));
            new_score = zeros(length(score)/num_time_idx,num_time_idx);
            for tt = 1:num_time_idx
                new_score(:,tt) = score(time_idx == tt);
            end

            if contains(fname,'not') == 1
                stats(1).time(1).freq(f).F.pt(pt_idx).not.data = spd.sim(f).F';
                stats(1).time(1).freq(f).score.pt(pt_idx).not.data = new_score;
            else
                stats(1).time(1).freq(f).F.pt(pt_idx).spike.data = spd.sim(f).F';
                stats(1).time(1).freq(f).score.pt(pt_idx).spike.data = new_score;
            end
        end

    end
    
    %% Add soz info
    if contains(fname,'not') == 0
        %is_soz_pt = zeros(length(ers.spike),1);
        is_soz_pt = zeros(length(spike.spike),1);
        for s = 1:length(spike.spike)
            if contains(met,'auto')
                biggest_dev = spike.spike(s).biggest_dev;
            else
                biggest_dev = manual_big.spike(s).dev_ch;
            end
            if ismember(biggest_dev,soz_chs)
                is_soz_pt(s) = 1;
            end
        end
        is_spike_soz(pt_idx).is_soz = is_soz_pt;
    end


end

    

%% Add spike power
%{
for n = 1:length(stats)
    for t = 1:length(stats(n).time)
        
        % confirm times align
        tw = stats(n).time(t).time_window;
        if tw ~= pre_spike(1).windows(t).which, error('what'); end
        
        for f = 1:length(stats(n).time(t).freq)
            
            
            %times = round((ers.index_windows(:,1)/ers.fs-3)*1e2)/(1e2);
            for p = 1:length(pre_spike)
                stats(n).time(t).freq(f).sd.pt(p).times = pre_spike(p).windows(t).cons_windows;
                %{
                if strcmp(wpr,'cons')
                    stats(n).time(t).freq(f).sd.pt(p).times = pre_spike(p).windows(t).cons_windows;
                else
                    stats(n).time(t).freq(f).sd.pt(p).times = pre_spike(p).windows(t).all_windows;
                end
                %}
                stats(n).time(t).freq(f).sd.pt(p).name = pre_spike(p).name;
                stats(n).time(t).freq(f).sd.pt(p).spike.data = pre_spike(p).windows(t).dev.spike;
                stats(n).time(t).freq(f).sd.pt(p).not.data = pre_spike(p).windows(t).dev.not;
            end
        end
    end
end
%}




end