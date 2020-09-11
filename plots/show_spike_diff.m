function show_spike_diff(windows)

%% Parameters
alpha = 0.05;

%% Get file locations, load spike times and pt structure
locations = spike_network_files;
main_folder = locations.main_folder;
results_folder = [main_folder,'results/'];
data_folder = [main_folder,'data/'];
eeg_folder = [main_folder,'results/eeg_data/'];
script_folder = locations.script_folder;
addpath(genpath(script_folder));
pt_file = [data_folder,'spike_structures/pt.mat'];
bct_folder = locations.BCT;
addpath(genpath(bct_folder));
out_folder = [results_folder,'plots/'];
sig_dev_folder = [results_folder,'signal_deviation/manual/'];
perm_folder = [results_folder,'perm_stats/'];
nbs_folder = [results_folder,'nbs_stats/'];
metric_folder = [results_folder,'metrics/manual/'];

freq_names = {'delta','theta','alpha','beta','low_gamma',...
    'high_gamma','ultra_high','broadband'};

%% Now get F statistics for network differences
network_count = 0;
n_freq_abs = 0;
listing = dir(perm_folder);
for l = 1:length(listing)
    name= listing(l).name;

    % Skip if . or ..
    if strcmp(name,'.') == 1 || strcmp(name,'..') == 1
        continue
    end

    % Skip if not a directory
    if listing(l).isdir == 0, continue; end

    network_count = network_count + 1;
    stats(network_count).name = name;

    network_folder = [perm_folder,name,'/'];

    % Loop through time scales
    time_listing = dir(network_folder);
    time_count = 0;

    for k = 1:length(time_listing)
        time_name= time_listing(k).name;
        time_window = str2num(time_name);

        % Skip if . or ..
        if strcmp(time_name,'.') == 1 || strcmp(time_name,'..') == 1
            continue
        end

        % Skip if not a directory
        if time_listing(k).isdir == 0, continue; end

        % Skip if not the time window we want
        if ismember(time_window,windows) == 0, continue; end

        time_count = time_count + 1;
        stats(network_count).time(time_count).name = time_name;
        stats(network_count).time(time_count).time_window = time_window;
        
        time_folder = [network_folder,time_name,'/'];

        pt_listing = dir([time_folder,'*.mat']);

        % load one to get nfreq
        sim = load([time_folder,pt_listing(1).name]);
        sim = sim.sim;
        nfreq = length(sim);
        if n_freq_abs < nfreq
            n_freq_abs = nfreq;
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
            sim = load([time_folder,fname]);
            sim = sim.sim;
            
            
            if isfield(sim(1),'index_windows') == 1 && isempty(sim(1).index_windows) == 0
                stats(network_count).time(time_count).index_windows = sim.index_windows;
                stats(network_count).time(time_count).fs = sim.fs;
            end


            for f = 1:nfreq
                stats(network_count).time(time_count).freq(f).name = sim(f).name;

                
                % spike vs not a spike
                if contains(fname,'not') == 0
                    stats(network_count).time(time_count).freq(f).F_all(pt_idx,:,1) = sim(f).F;
                else
                    stats(network_count).time(time_count).freq(f).F_all(pt_idx,:,2) = sim(f).F;
                end
            end


        end

    end

end

%% Show F stats for spike vs not a spike at the spike time
n = 2; % correlation
F = squeeze(stats(n).time(1).freq.F_all(:,size(stats(n).time(1).freq.F_all,2),:));
figure
plot(1,F(:,1),'ko','MarkerSize',18)
hold on
plot(2,F(:,2),'ko','MarkerSize',18)
set(gca,'xlim',[0.5 2.5])

% Do sign rank test
p = signrank(F(:,1),F(:,2));
if p < 0.001
    p_text = sprintf('p < 0.001');
else
    p_text = sprintf('p = %1.3f',p);
end

plot([1 2],[max(max(F)) max(max(F))] + 2,'k','linewidth',2)
text(1.5,max(max(F))+4,p_text,'fontsize',20,'horizontalalignment','center')
set(gca,'fontsize',20)
xticks([1 2])
xticklabels({'Spike','Not spike'})
ylabel('pseudo-F statistic')
ylim([0 max(max(F))+6])

print(gcf,[out_folder,'spike_vs_no_F'],'-depsc')

%% Now get ns and ge for spike vs no spike
network_count = 0;
n_freq_abs = 0;
listing = dir(metric_folder);
clear stats
for l = 1:length(listing)
    name= listing(l).name;

    % Skip if . or ..
    if strcmp(name,'.') == 1 || strcmp(name,'..') == 1
        continue
    end

    % Skip if not a directory
    if listing(l).isdir == 0, continue; end

    network_count = network_count + 1;
    stats(network_count).name = name;

    network_folder = [metric_folder,name,'/'];

    % Loop through time scales
    time_listing = dir(network_folder);
    time_count = 0;

    for k = 1:length(time_listing)
        time_name= time_listing(k).name;
        time_window = str2num(time_name);

        % Skip if . or ..
        if strcmp(time_name,'.') == 1 || strcmp(time_name,'..') == 1
            continue
        end

        % Skip if not a directory
        if time_listing(k).isdir == 0, continue; end

        % Skip if not the time window we want
        if ismember(time_window,windows) == 0, continue; end

        time_count = time_count + 1;
        stats(network_count).time(time_count).name = time_name;
        stats(network_count).time(time_count).time_window = time_window;
        
        time_folder = [network_folder,time_name,'/'];

        pt_listing = dir([time_folder,'*.mat']);

        % load one to get nfreq
        sim = load([time_folder,pt_listing(1).name]);
        sim = sim.metrics;
        nfreq = length(sim.freq);
        if n_freq_abs < nfreq
            n_freq_abs = nfreq;
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
            sim = load([time_folder,fname]);
            sim = sim.metrics;
            
            
            if isfield(sim(1),'index_windows') == 1 && isempty(sim(1).index_windows) == 0
                stats(network_count).time(time_count).index_windows = sim.index_windows;
                stats(network_count).time(time_count).fs = sim.fs;
            end


            for f = 1:nfreq
                stats(network_count).time(time_count).freq(f).name = sim.freq(f).name;


                stats(network_count).time(time_count).freq(f).ns_all(pt_idx,:) =...
                    mean(sim.freq(f).ns.data,1);
                
                stats(network_count).time(time_count).freq(f).ge_all(pt_idx,:) =...
                    mean(sim.freq(f).ge.data,1);

            end


        end

    end

end

%% Show GE for spike vs not a spike at the spike time
n = 2; % correlation
ge = stats(n).time(1).freq.ge_all(:,[1,size(stats(n).time(1).freq.ge_all,2)]);
figure
plot(1,ge(:,1),'ko','MarkerSize',18)
hold on
plot(2,ge(:,2),'ko','MarkerSize',18)
set(gca,'xlim',[0.5 2.5])

% Do sign rank test
p = signrank(ge(:,1),ge(:,2));
if p < 0.001
    p_text = sprintf('p < 0.001');
else
    p_text = sprintf('p = %1.3f',p);
end

plot([1 2],[max(max(ge)) max(max(ge))] + 2,'k','linewidth',2)
text(1.5,max(max(ge))+4,p_text,'fontsize',20,'horizontalalignment','center')
set(gca,'fontsize',20)
xticks([1 2])
xticklabels({'Pre-spike','Spike'})
ylabel('Global efficiency')
ylim([0 max(max(ge))+6])

%print(gcf,[out_folder,'spike_vs_no_F'],'-depsc')

end