function stats = get_specified_metrics(windows,pre_spike,met)

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

%% Now get F statistics for network differences
time_name = sprintf('%1.1f/',windows);

if strcmp(met,'ers')
    main_folder = ers_folder;
elseif strcmp(met,'F')
    main_folder =  sp_diff_folder;
elseif strcmp(met,'ns_inv') || strcmp(met,'ge') || strcmp(met,'ns_all') || ...
                strcmp(met,'ns_avg') || strcmp(met,'trans')
    main_folder = ns_folder;
end

time_folder = [main_folder,time_name];

pt_listing = dir([time_folder,'*.mat']);


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
    load([time_folder,fname]);


    % Get number of frequencies
    if strcmp(met,'ers')
        nfreq = size(ers.freq_bands,1);
        for f = 1:nfreq
            stats(1).time(1).freq(f).name = ers.freq_names{f};

        end
    elseif strcmp(met,'F')
        nfreq = length(sim);
        for f = 1:nfreq
            stats(1).time(1).freq(f).name = sim.name;
        end
    elseif strcmp(met,'ns_inv') || strcmp(met,'ge') || strcmp(met,'ns_all') || ...
                strcmp(met,'ns_avg') || strcmp(met,'trans')
        nfreq = length(metrics.freq);
        for f = 1:nfreq
            stats(1).time(1).freq(f).name = metrics.freq(f).name;
        end
    else

        nfreq = length(sim);
        for f = 1:nfreq
            stats(1).time(1).freq(f).name = sim.freq(f).name;
        end
    end
    stats(1).time(1).time_window = windows;

    for f = 1:nfreq
        if strcmp(met,'ns_inv') || strcmp(met,'ge') || strcmp(met,'ns_all') || ...
                strcmp(met,'ns_avg') || strcmp(met,'trans')

            sim = metrics;

            ns = sim.freq(f).ns.data;
            ge = sim.freq(f).ge.data;
            ns_all = sim.freq(f).ns_all.data;
            trans = sim.freq(f).trans.data;
            involved = logical(sim.freq(f).involved);

            % Avg ns_all across channels
            ns_avg = squeeze(mean(ns_all,3));

            % Average ns_inv across involved channels
            ns_inv = zeros(size(ns_avg));
            for s = 1:size(ns_all,1)
                curr_sp = squeeze(ns_all(s,:,:));
                curr_inv = involved(s,:);
                curr_ns_inv = squeeze(mean(curr_sp(:,curr_inv),2));
                ns_inv(s,:) = curr_ns_inv;
            end

            times = round((sim.index_windows(:,1)/sim.fs-3)*1e2)/(1e2);
            
            stats(1).time(1).freq(f).ge.pt(pt_idx).name = pt_name;
            stats(1).time(1).freq(f).ns_avg.pt(pt_idx).name = pt_name;
            stats(1).time(1).freq(f).ns_big.pt(pt_idx).name = pt_name;
            stats(1).time(1).freq(f).ns_inv.pt(pt_idx).name = pt_name;
            stats(1).time(1).freq(f).trans.pt(pt_idx).name = pt_name;

            % spike vs not a spike
            if contains(fname,'not') == 1
                stats(1).time(1).freq(f).ge.pt(pt_idx).not.data(:,:) = ge;
                stats(1).time(1).freq(f).ns_avg.pt(pt_idx).not.data(:,:) = ns_avg;
                stats(1).time(1).freq(f).ns_big.pt(pt_idx).not.data(:,:) = ns;
                stats(1).time(1).freq(f).ns_inv.pt(pt_idx).not.data(:,:) = ns_avg; % just average for not spike
                stats(1).time(1).freq(f).trans.pt(pt_idx).not.data(:,:) = trans;
            else
                stats(1).time(1).freq(f).ge.pt(pt_idx).spike.data(:,:) = ge;
                stats(1).time(1).freq(f).ns_avg.pt(pt_idx).spike.data(:,:) = ns_avg;
                stats(1).time(1).freq(f).ns_big.pt(pt_idx).spike.data(:,:) = ns;
                stats(1).time(1).freq(f).ns_inv.pt(pt_idx).spike.data(:,:) = ns_inv; 
                stats(1).time(1).freq(f).trans.pt(pt_idx).spike.data(:,:) = trans;
            end

            stats(1).time(1).freq(f).ns_big.name = 'Node strength (spike channel)';
            stats(1).time(1).freq(f).ns_inv.name = 'Node strength (involved channels)';
            stats(1).time(1).freq(f).ns_avg.name = 'Node strength (average)';
            stats(1).time(1).freq(f).ge.name = 'Global efficiency';
            stats(1).time(1).freq(f).trans.name = 'Transitivity';

            stats(1).time(1).freq(f).ns_big.pt(pt_idx).times = times;
            stats(1).time(1).freq(f).ns_inv.pt(pt_idx).times = times;
            stats(1).time(1).freq(f).ns_avg.pt(pt_idx).times = times;
            stats(1).time(1).freq(f).ge.pt(pt_idx).times = times;
            stats(1).time(1).freq(f).trans.pt(pt_idx).times = times;
        end

        % Add ERS stuff
        if strcmp(met,'ers')

            stats(1).time(1).freq(f).ers.pt(pt_idx).index_windows = ers.index_windows;
            times = round((ers.index_windows(:,1)/ers.fs-3)*1e2)/(1e2);
            stats(1).time(1).freq(f).ers.pt(pt_idx).times = times;
            stats(1).time(1).freq(f).ers.pt(pt_idx).name = ers.name;
            
            %
            all_ers(f).data = zeros(length(ers.spike),size(ers.spike(1).ers,1));
            for s = 1:length(ers.spike)
                all_ers(f).data(s,:) = ers.spike(s).ers(:,f);
            end
           

            if contains(fname,'not') == 1
                stats(1).time(1).freq(f).ers.pt(pt_idx).not.data = all_ers(f).data;
            else
                stats(1).time(1).freq(f).ers.pt(pt_idx).spike.data = all_ers(f).data;
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


end

    

%% Add spike power
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


end