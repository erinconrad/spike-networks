function ers_ch_comparison

%% Parameters
alpha = 0.05;
which_freq = 4;
rm_rise = 1; 
met = 'ers_all';
windows = [0.1];
method = 'ttestp'; % ttestp is default
which_pt = 1;
which_pre_rise = 0; % 2 is default
comp_points = 3;  %3 is default
% 0 = absolute, 1 = z score, 2 = relative change from first one, 3 = like z
% score but subtracting first one

if which_pre_rise == 0
    wpr = 'manual_before_rise';
elseif which_pre_rise == 1
    wpr = 'before_rise';
    error('I do not do this anymore')
else
    wpr = 'cons';
end

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
adj_folder = [results_folder,'adj_mat/manual/'];
spike_rise_folder = [results_folder,'spike_rise/'];
pre_spike_folder = [results_folder,'/pre_spike/'];

if exist(out_folder,'dir') == 0
    mkdir(out_folder);
end


if exist(pre_spike_folder,'dir') == 0
    mkdir(pre_spike_folder);
end


%% Find windows for each spike in which no early spike rise
pre_spike = find_pre_spike_windows(windows);

%% Get signal power deviation
sig_dev = get_sd(alpha,0,0);
% convert this to be similar to pre_spike
pre_spike = convert_sd(sig_dev,windows,pre_spike);

%% Get ERS
ers_subfolder = [ers_folder,sprintf('%1.1f',windows),'/'];
listing = dir([ers_subfolder,'*ers.mat']);
all_names = {};
for l = 1:length(listing)
    fname= listing(l).name;
    
    pt_name = strsplit(fname,'_');
    pt_name = pt_name{1};


    [a,b] = ismember(pt_name,all_names);
    if a == 1
        pt_idx = b;
    else
        all_names = [all_names;pt_name];
        pt_idx = length(all_names);
    end
    
    % Load file
    ers_file = [ers_subfolder,fname];
    ers = load(ers_file);
    ers = ers.ers;
    nfreq = size(ers.freq_bands,1);
    
    for f = 1:nfreq
        freq(f).ers.pt(pt_idx).index_windows = ers.index_windows;
        times = round((ers.index_windows(:,1)/ers.fs-3)*1e2)/(1e2);
        freq(f).ers.pt(pt_idx).times = times;
        freq(f).ers.pt(pt_idx).name = ers.name;

        % Concatenate spike data
        all_ers(f).data = zeros(length(ers.spike),size(ers.spike(1).ers_all,1),size(ers.spike(1).ers_all,3));
        for s = 1:length(ers.spike)       
            curr_ers = squeeze(ers.spike(s).ers_all(:,f,:));
            all_ers(f).data(s,:,:) = curr_ers;
            
        end
        
        if contains(fname,'not') == 1
            freq(f).ers.pt(pt_idx).not.data = all_ers(f).data;
        else
            freq(f).ers.pt(pt_idx).spike.data = all_ers(f).data;
        end
        
    end
    
end

%% Get the last non-removed time
for f = 1:length(freq)
    curr_met = freq(f).ers;
    
    if strcmp(wpr,'cons')
        ps_windows = pre_spike(1).windows(1).cons_windows;
    else
        ps_windows = pre_spike(1).windows(1).all_windows;
    end
    shift = find(ps_windows == curr_met.pt(1).times(1));

    all_t = nan(length(curr_met.pt),length(curr_met.pt(1).times));
    
    % Get early rise times
    if strcmp(wpr,'cons')
        % Loop over patients
        for p = 1:length(curr_met.pt)
            spike_dev = pre_spike(p).windows(1).dev.spike(:,shift:end);

            % Paired ttest comparing first time window
            % to subsequent time windows
            for tt = 2:size(spike_dev,2)
                [~,~,~,tstats] = ttest(spike_dev(:,1),spike_dev(:,tt));
                all_t(p,tt) = tstats.tstat;

            end

        end

        % Do one sample test of tstats
        before_rise_windows = ones(size(all_t,2),1);
        for tt = 2:size(all_t,2)
            [~,pval] = ttest(all_t(:,tt));
            if pval < alpha
                before_rise_windows(tt:end) = 0;
            end
        end

    end
    
    % Loop over patients
    for p = 1:length(curr_met.pt)

        if p > length(curr_met.pt), continue; end

        curr_pt = curr_met.pt(p);

        if ~isfield(curr_pt,'spike') || ~isfield(curr_pt,'not')
            continue;
        end

        % Skip it if I don't have pre-spike info for it
        if p > length(pre_spike), continue; end

        % check that I have the same pt in pre_spike for
        % this index
        if ~strcmp(pre_spike(p).name,curr_pt.name)
            error('Names do not align');
        end

        if strcmp(wpr,'cons')
            before_rise = repmat(before_rise_windows',...
                size(curr_pt.spike.data,1),1);
        else
            before_rise = pre_spike(p).windows(1).(wpr);
        end
        
        if size(before_rise,1) > size(curr_pt.spike.data,1)
            before_rise = before_rise(1:size(curr_pt.spike.data,1),:);
        end

        if size(before_rise,2) > size(curr_pt.spike.data,2)
            before_rise = before_rise(:,size(before_rise,2) - size(curr_pt.spike.data,2) +1: end);
        end
        
        % Get the last time window before rise
        last_before_rise = zeros(size(before_rise,1),1);
        for l = 1:size(before_rise,1)
            curr_spike_rise = before_rise(l,:);
            last = find(curr_spike_rise == 1,1,'last');
            last_before_rise(l) = last;
        end
        
        % Get the ERS data from this last time
        last_sp = zeros(size(curr_pt.spike.data,1),size(curr_pt.spike.data,3));
        last_not = zeros(size(curr_pt.not.data,1),size(curr_pt.not.data,3));
        
        for s = 1:size(curr_pt.spike.data,1)
            last_sp(s,:) = squeeze(curr_pt.spike.data(s,last_before_rise(s),:));
        end
        for s = 1:size(curr_pt.not.data,1)
            last_not(s,:) = squeeze(curr_pt.not.data(s,last_before_rise(s),:));
        end
        freq(f).ers.pt(pt_idx).spike.last = last_sp;
        freq(f).ers.pt(pt_idx).not.last = last_not;
        
        % Relative power in biggest dev ch
        freq(f).ers.pt(pt_idx).spike.big_dev = 
        
    end
end


end