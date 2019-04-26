function calc_network_metrics(whichPts)

%% Parameters
% 1 = alpha/theta; 2 = beta, 3 = low gamma, 4 = high gamma, 5 = ultra high, 6 = broadband
freq_text = {'alpha/theta','beta','low\ngamma','high\ngamma','ultra high\ngamma','broadband'};
n_f = length(freq_text);
n_times = 11;
spike_window_times = [-0.2 0.8];

%% Get file locations, load spike times and pt structure
locations = spike_network_files;
main_folder = locations.main_folder;
results_folder = [main_folder,'results/'];
data_folder = [main_folder,'data/'];
script_folder = locations.script_folder;
addpath(genpath(script_folder));
spike_times_file = [data_folder,'spike_times/times.mat'];
pt_file = [data_folder,'spike_structures/pt.mat'];
bct_folder = locations.BCT;
addpath(genpath(bct_folder));
plot_folder = [results_folder,'plots/network_dev/'];
if exist(plot_folder,'dir') == 0
    mkdir(plot_folder);
end

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
    adj_folder = [results_folder,name,'/adj/'];
    fs = pt(whichPt).fs;
    
    %% Load adjacency matrices and calculate metrics
    listing = dir([adj_folder,'adj*.mat']);
    if length(listing) == 0
        fprintf('No adjacency matrices for %s, skipping...\n',name);
        continue
    end
    
    % Number of spikes expected
    n_spikes = sum(~isnan(times(whichPt).spike_times));
    
    % Prep network matrices
    ge = nan(n_f,n_spikes,n_times);
    ns_seq = nan(n_f,n_spikes,n_times);
    sync = nan(n_f,n_spikes,n_times);
    ec_seq = nan(n_f,n_spikes,n_times);
    dev = nan(n_spikes,fs*(n_times+1));
    bin_dev = nan(n_spikes,n_times);
    
    % Initialize spike count
    s_count = 0;
    
    for f = 1:length(listing)
        fname = listing(f).name;
        meta = load([adj_folder,fname]); %outputs a struct named meta
        
        % Load the adjacency matrix file
        meta = meta.meta;
        
        % Load the spike file
        spike = load([pt_folder,sprintf('spikes_%d.mat',f)]);
        spike = spike.spike;
        
        % Loop through spikes
        for s = 1:length(meta.spike)

            s_count = s_count + 1;
            
            fprintf('Doing spike %d of %d\n',s_count,length(listing)*100);
            
            
            seq_chs = meta.spike(s).is_seq_ch; % binary array
            sp_ch = meta.spike(s).is_sp_ch; % binary array
            
            %% Get spike signal
            data = load([pt_folder,'basic_info.mat']); % returns a structure called data
            data = data.data;
            
            ch_labels = data.chLabels(:,1);
            ignore = zeros(length(ch_labels),1);

            for i = 1:length(ch_labels)
                ch_labels{i} = ieeg_ch_parser(ch_labels{i});
                for j = 1:length(pt(whichPt).ignore.names)
                    if strcmp(pt(whichPt).ignore.names(j),ch_labels{i}) == 1
                        ignore(i) = 1;
                    end
                end
            end
            
            % Get the values
            values = spike(s).values(:,~ignore);
            values = values(:,sp_ch);
            
            % Divide it up into windows
            peak = round(size(values,1)/2);
            index_windows = zeros(n_times,2); %6 is the spike window
            tick_window = 1*fs;
            spike_window = peak + spike_window_times*fs;
            
            for i = 1:size(index_windows,1)
                
                % if we're at the spike window (6), then the index window
                % is the spike window. If i is 1, then it's 5 seconds
                % before. If i is 11, it's 5 seconds after
                index_windows(i,:) = spike_window + tick_window*(i-6);
            end
            
            % readjust size
            if length(values) < size(dev,2)
                fprintf('Padding values by %d\n',size(dev,2)-length(values));
                values = [values;nan(size(dev,2)-length(values),1)];
                
            elseif length(values) > size(dev,2)
                fprintf('Shortening values by %d\n',size(dev,2)-length(values));
                values = values(1:end-(length(values)-size(dev,2)));
                
            end
            
            dev_t = sqrt((values-nanmedian(values)).^2);
            dev(s_count,:) = dev_t;
            
            for i = 1:size(index_windows,1)
                new_dev = sqrt((values-nanmedian(values)).^2);
                bin_dev(s_count,i) = nanmean(new_dev(round(index_windows(i,1)):...
                    round(index_windows(i,2))));
            end
            
            for which_freq = 1:length(meta.spike(s).adj)
                adj_all_t= meta.spike(s).adj(which_freq).adj;

                % Loop through times
                for tt = 1:size(adj_all_t,1)

                    % Get adj matrix of interest
                    adj = squeeze(adj_all_t(tt,:,:));

                    %% Calculate metrics
                    ge(which_freq,s_count,tt) = efficiency_wei(adj,0); % global efficiency
                    sync(which_freq,s_count,tt) = synchronizability_sp(adj);

                    % node strength of all ch in seq
                    ns_temp = strengths_und(adj); 
                    ns_seq(which_freq,s_count,tt) = sum(ns_temp(seq_chs));

                    % eigenvector centrality of all ch in seq
                    ec_temp = eigenvector_centrality_und(adj);
                    ec_seq(which_freq,s_count,tt) = sum(ec_temp(seq_chs));

                end
                
            end
            
        end
        
    end
    
    %% Aggregate metrics
    avg_ns_seq = nanmean(ns_seq-median(ns_seq,3),2);
    avg_ge = nanmean(ge - median(ge,3),2);
    avg_sync = nanmean(sync - median(sync,3),2);
    avg_ec_seq = nanmean(ec_seq - median(ec_seq,3),2);
    
    plot_thing(1,:,:) = avg_ns_seq;
    plot_thing(2,:,:) = avg_ec_seq;
    plot_thing(3,:,:) = avg_ge;
    plot_thing(4,:,:) = avg_sync;
    
    plot_title{1} = 'Node strength\nof spike sequence chs';
    plot_title{2} = 'Eigenvector centrality\nof spike sequence chs';
    plot_title{3} = 'Global efficiency';
    plot_title{4} = 'Synchronizability';
    
    %% Get avg deviation of signal
    avg_dev = nanmean(dev,1);
    avg_bin_dev = nanmean(bin_dev,1);
    
    %% Plot aggregated metrics
    figure
    set(gcf,'position',[26 0 1242 900])
    [ha, pos] = tight_subplot(n_f, 4, [0.04 0.04], [0.08 0.08], [0.05 0.01]);
    for f = 1:n_f
        for i = 1:size(plot_thing,1)
            axes(ha((f-1)*4+i))
            plot(squeeze(plot_thing(i,f,:)),'ks-','linewidth',2)
            if f == 1
                title(sprintf(plot_title{i}));
            end
            if f == n_f
               xlabel('Time (s)') 
            end
            yticklabels([])
            
            if i == 1
                ylabel(sprintf(freq_text{f}));
            end
            set(gca,'fontsize',20)
        end
    end
    filename = [name,'_network_dev'];
    print([plot_folder,filename],'-depsc');
    
    
    figure
    set(gcf,'position',[26 0 1100 400])
    [ha2, pos] = tight_subplot(1, 2, [0.01 0.02], [0.14 0.08], [0.01 0.01]);
    axes(ha2(1))
    plot((1:length(avg_dev))/fs,avg_dev,'k-')
    hold on
    for tt = 1:size(index_windows,1)
        plot([index_windows(tt,1) index_windows(tt,1)]/fs,get(gca,'ylim'),'k--')
        plot([index_windows(tt,2) index_windows(tt,2)]/fs,get(gca,'ylim'),'k--')
        text((index_windows(tt,1)+index_windows(tt,2))/2/fs-0.25,max(avg_dev),...
            sprintf('%d',tt),'fontsize',20);
    end
    yticklabels([])
    xlabel('Time (s)')
    title('Average signal deviation from baseline')
    set(gca,'fontsize',20)
    
    axes(ha2(2))
    plot(avg_bin_dev,'ks-','linewidth',2)
    title('Average binned signal deviation from baseline')
    yticklabels([])
    xlabel('Time (s)')
    set(gca,'fontsize',20)
    
    filename = [name,'_signal_dev'];
    print([plot_folder,filename],'-depsc');
    
   
    
end


end