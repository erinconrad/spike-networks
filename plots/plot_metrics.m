function plot_metrics(t_window)

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
metrics_folder = [results_folder,'metrics/manual/'];


if exist(out_folder,'dir') == 0
    mkdir(out_folder);
end

% Load spike file for one patient to get the surround time
spike = load([eeg_folder,'HUP074_eeg.mat']);
spike = spike.spike;
surround_time = spike(1).surround_time;

% Loop through network types
listing = dir(metrics_folder);
network_count = 0;
n_freq_abs = 0;
max_F = 0;
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
    
    network_folder = [metrics_folder,name,'/'];
    
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
        
        % skip if not the one i want
        if t_window ~= time_window
            continue
        end
        
        time_count = time_count + 1;
        stats(network_count).time(time_count).name = time_name;
        stats(network_count).time(time_count).time_window = time_window;
        time_folder = [network_folder,time_name,'/'];
        
        pt_listing = dir([time_folder,'*.mat']);
        
        % load one to get nfreq
        metrics = load([time_folder,pt_listing(1).name]);
        metrics = metrics.metrics;
        nfreq = length(metrics.freq);
        if n_freq_abs < nfreq
            n_freq_abs = nfreq;
        end
        
        for f = 1:nfreq
            stats(network_count).time(time_count).freq(f).ge.data = ...
                zeros(length(pt_listing),surround_time*2/time_window);

            stats(network_count).time(time_count).freq(f).p_all = ...
                zeros(length(pt_listing),surround_time*2/time_window,2);
        end
        
        % loop through pts
        for i = 1:length(pt_listing)
            
            pt_name = pt_listing(i).name;
            pt_name_pt = strsplit(pt_name,'_');
            pt_name_pt = pt_name_pt{1};
            
            % load pt file
            metrics = load([time_folder,pt_name]);
            metrics = metrics.metrics;
            involved = metrics.involved;
            
            for f = 1:nfreq
                stats(network_count).time(time_count).freq(f).name = metrics.freq(f).name;
                stats(network_count).time(time_count).freq(f).ge.name = metrics.freq(f).ge.name;
                stats(network_count).time(time_count).freq(f).ge.data(i,:) = metrics.freq(f).ge.data;
                stats(network_count).time(time_count).freq(f).ns.name = metrics.freq(f).ns.name;
                stats(network_count).time(time_count).freq(f).bc.name = metrics.freq(f).bc.name;
                
                % Get ns, bc for involved and uninvolved chs (1 = involved,
                % 2 = uninvolved)
                stats(network_count).time(time_count).freq(f).ns.data(i,:,1) = ...
                    mean(metrics.freq(f).ns.data(:,involved),2);
                stats(network_count).time(time_count).freq(f).ns.data(i,:,2) = ...
                    mean(metrics.freq(f).ns.data(:,~involved),2);
                
                stats(network_count).time(time_count).freq(f).bc.data(i,:,1) = ...
                    mean(metrics.freq(f).bc.data(:,involved),2);
                stats(network_count).time(time_count).freq(f).bc.data(i,:,2) = ...
                    mean(metrics.freq(f).bc.data(:,~involved),2);

            
            end
            
        end
        
    end
 
end

%% Plot
colors = [0, 0.4470, 0.7410;0.8500, 0.3250, 0.0980;0.9290, 0.6940, 0.1250];

% Initialize figure
%{
nfreq (coherence) + 1 (simple) columns and 2 (2 time scales) rows
%}
figure
set(gcf,'position',[100 100 1300 500])
[ha, pos] = tight_subplot(3, n_freq_abs+1, [0.01 0.01], [0.1 0.05], [0.05 0.01]);

for n = 1:network_count

    net_name = stats(n).name;
    
    for t = 1:time_count
        
        % change times for x axis
        nchunks = size(stats(n).time(t).freq(1).ge.data,2);
        times = realign_times(nchunks,surround_time);
        
        nfreq = length(stats(n).time(t).freq);
        for f = 1:nfreq
            
            % Get appropriate subplot
            if strcmp(net_name,'coherence') == 1
                column_add = 1;
            else
                column_add = 0;
            end
            
          
            % get mean and std across pts for each of the network metrics
            ge_mean = nanmean(stats(n).time(t).freq(f).ge.data,1);
            ns_in_mean = nanmean(stats(n).time(t).freq(f).ns.data(:,:,1),1);
            bc_in_mean = nanmean(stats(n).time(t).freq(f).ns.data(:,:,1),1);
            
            ns_out_mean = nanmean(stats(n).time(t).freq(f).ns.data(:,:,2),1);
            bc_out_mean = nanmean(stats(n).time(t).freq(f).ns.data(:,:,2),1);
            
            ge_std = nanstd(stats(n).time(t).freq(f).ge.data,1);
            ns_in_std = nanstd(stats(n).time(t).freq(f).ns.data(:,:,1),1);
            bc_in_std = nanstd(stats(n).time(t).freq(f).ns.data(:,:,1),1);
            
            ns_out_std = nanstd(stats(n).time(t).freq(f).ns.data(:,:,2),1);
            bc_out_std = nanstd(stats(n).time(t).freq(f).ns.data(:,:,2),1);
            
            max_val = max([ge_mean,ns_in_mean,bc_in_mean,ns_out_mean,bc_out_mean]);
            
            % plot means with stds as error bars
            sp = f + column_add;
            axes(ha(sp));
            hold on
            gep = errorbar(times,ge_mean,ge_std,'color',colors(1,:));
            hold on
            
            if strcmp(net_name,'coherence') == 1
                title(sprintf('%s',...
                    strrep(stats(n).time(t).freq(f).name,'_',' ')))
            elseif strcmp(net_name,'simple') == 1
                title('correlation')
            end
            
            sp = (n_freq_abs+1) + f + column_add;
            axes(ha(sp));
            hold on
            nsip = errorbar(times,ns_in_mean,ns_in_std,'color',colors(2,:));
            hold on
            %nsop = errorbar(times,ns_out_mean,ns_out_std,'--','color',colors(2,:));
            
            
            sp = (n_freq_abs+1)*2 + f + column_add;
            axes(ha(sp));
            hold on
            bcip = errorbar(times,bc_in_mean,bc_in_std,'color',colors(3,:));
            hold on
            %bcop = errorbar(times,bc_out_mean,bc_out_std,'--','color',colors(3,:));
            
            if f == 4
                xlabel('Time relative to spike peak (s)')
            end 
            
            
            % Determine significance by paired t-test comparing first
            % second to current second
            ge_p = nan(length(ge_mean),1);
            ns_p = nan(length(ge_mean),1);
            bc_p = nan(length(ge_mean),1);
            for tt = 2:length(ge_mean)
                
                [~,ge_p(tt)] = ttest(stats(n).time(t).freq(f).ge.data(:,1,1),...
                    stats(n).time(t).freq(f).ge.data(:,tt,1));
                [~,ns_p(tt)] = ttest(stats(n).time(t).freq(f).ns.data(:,1,1),...
                    stats(n).time(t).freq(f).ns.data(:,tt,1));
                [~,bc_p(tt)] = ttest(stats(n).time(t).freq(f).bc.data(:,1,1),...
                    stats(n).time(t).freq(f).bc.data(:,tt,1));
            end
            
            % Display asterisks if significant p-values
            for tt = 2:length(ge_p)
                sp = f + column_add;
                axes(ha(sp));
                text_out = get_asterisks(ge_p(tt),nchunks*(n_freq_abs+1));
                text(times(tt),ge_mean(tt)+ge_std(tt)+0.5,...
                    sprintf('%s',text_out),'fontsize',20,'horizontalalignment','center')
               
                sp = (n_freq_abs+1) + f + column_add;
                axes(ha(sp));
                text_out = get_asterisks(ns_p(tt),nchunks*(n_freq_abs+1));
                text(times(tt),ns_in_mean(tt)+ns_in_std(tt)+0.5,...
                    sprintf('%s',text_out),'fontsize',20,'horizontalalignment','center')
                
                sp = (n_freq_abs+1)*2 + f + column_add;
                axes(ha(sp));
                text_out = get_asterisks(bc_p(tt),nchunks*(n_freq_abs+1));
                text(times(tt),bc_in_mean(tt)+bc_in_std(tt)+0.5,...
                    sprintf('%s',text_out),'fontsize',20,'horizontalalignment','center')
            end
           
           
            
        end
        
        
        
    end
end

%{
for sp = 1:length(ha)
    axes(ha(sp))
    % formatting
    ylim([0 max_F+2]);
    if mod(sp,n_freq_abs+1) == 1
        ylabel(sprintf('%s s time window',stats(1).time(floor((sp/n_freq_abs+1))).name));
    else
        yticklabels([])
    end
    
    set(gca,'fontsize',15)
end
%}

print(gcf,[out_folder,'metric_change'],'-depsc');


end