function just_spike_analysis(whichPts)

freq_text = {'alpha/theta','beta','low\ngamma','high\ngamma','ultra high\ngamma','broadband'};
%freq_text = {'alpha/theta'};
n_f = length(freq_text);
n_times = 3;
n_seconds = 14;

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
plot_folder = [results_folder,'plots/just_spike/'];

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
    adj_folder = [results_folder,name,'/just_spike/'];
    fs = pt(whichPt).fs;
    
    stats_folder = [pt_folder,'stats/'];
    if exist(stats_folder,'dir') == 0
        mkdir(stats_folder);
    end
    
    %% Load adjacency matrices and calculate metrics
    listing = dir([adj_folder,'adj*.mat']);
    if length(listing) == 0
        fprintf('No adjacency matrices for %s, skipping...\n',name);
        continue
    end
    
    % Number of spikes expected
    n_spikes = sum(~isnan(times(whichPt).spike_times));
    
    % Load the first adjacency matrix to get the size in order to
    % initialize one of the arrays
    meta = load([adj_folder,listing(1).name]); 
    meta = meta.meta;
    nchs = length(meta.spike(1).is_seq_ch);
    
    % Prep network matrices
    ge = nan(n_f,n_spikes,n_times);
    ns_seq = nan(n_f,n_spikes,n_times);
    ns_not_seq = nan(n_f,n_spikes,n_times);
    sync = nan(n_f,n_spikes,n_times);
    ec_seq = nan(n_f,n_spikes,n_times);
    ec_not_seq = nan(n_f,n_spikes,n_times);
    dev = nan(n_spikes,fs*n_seconds);
    bin_dev = nan(n_spikes,n_times);
    is_seq = nan(n_spikes,nchs);
    avg_adj = zeros(n_f,n_times,nchs,nchs);
    
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
            is_seq(s,:) = seq_chs;
            
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
            index_windows = meta.spike(s).windows;
            
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
                bin_dev(s_count,i) = nanmean(dev_t(round(index_windows(i,1)):...
                    round(index_windows(i,2))));
            end
            
            for which_freq = 1:length(meta.spike(s).adj)
                adj_all_t= meta.spike(s).adj(which_freq).adj;

                % Loop through times
                for tt = 1:size(adj_all_t,1)

                    % Get adj matrix of interest
                    adj = squeeze(adj_all_t(tt,:,:));
                    
                    
                    
                    
                    if sum(sum(isnan(adj))) > 0
                        continue
                    end
                    
                    avg_adj(which_freq,tt,:,:) = squeeze(avg_adj(which_freq,tt,:,:)) + ...
                        adj;
                    
                    %% Calculate metrics
                    ge(which_freq,s_count,tt) = efficiency_wei(adj,0); % global efficiency
                    sync(which_freq,s_count,tt) = synchronizability_sp(adj);

                    % node strength of all ch in seq
                    ns_temp = strengths_und(adj); 
                    ns_seq(which_freq,s_count,tt) = mean(ns_temp(seq_chs));
                    ns_not_seq(which_freq,s_count,tt) = mean(ns_temp(~seq_chs));

                    % eigenvector centrality of all ch in seq
                    ec_temp = eigenvector_centrality_und(adj);
                    ec_seq(which_freq,s_count,tt) = mean(ec_temp(seq_chs));
                    
                    ec_not_seq(which_freq,s_count,tt) = mean(ec_temp(~seq_chs));

                end
                
            end
            
        end
        
    end
    
    %% Aggregate metrics
    
    avg_adj = avg_adj/s_count;
    
    avg_ns = nanmean(ns_seq-median(ns_seq,3),2);
    avg_ge = nanmean(ge - median(ge,3),2);
    avg_sync = nanmean(sync - median(sync,3),2);
    avg_ec = nanmean(ec_seq - median(ec_seq,3),2);
    %}
    
    %% Save stuff to structure
    out = [];
    out.signal.dev = dev;
    out.signal.bin_dev = bin_dev;
    
    out.network.ns = ns_seq;
    out.network.ec = ec_seq;
    out.network.ns_notseq = ns_not_seq;
    out.network.ec_notseq = ec_not_seq;
    out.network.ge = ge;
    out.network.sync = sync;
    
    avg_dev = nanmean(dev,1);
    avg_bin_dev = nanmean(bin_dev,1);
    
    save([stats_folder,'spike_stats.mat'],'out');
    
    plot_thing(1,:,:) = avg_ns;
    plot_thing(2,:,:) = avg_ec;
    plot_thing(3,:,:) = avg_ge;
    plot_thing(4,:,:) = avg_sync;
    
    plot_title{1} = 'Node strength\nof spike sequence chs';
    plot_title{2} = 'Eigenvector centrality\nof spike sequence chs';
    plot_title{3} = 'Global efficiency';
    plot_title{4} = 'Synchronizability';
    
    
    figure
    set(gcf,'position',[240 1 953 800])
    [ha, pos] = tight_subplot(n_f, n_times, [0.01 0.02], [0 0.04], [0.01 0.01]);
    for f = 1:n_f
        for i = 1:n_times
            axes(ha((f-1)*n_times+i))
            imagesc(squeeze(avg_adj(f,i,:,:)))
            xticklabels([])
            yticklabels([])
            colorbar
            
            if f == 1
                title(sprintf('Window %s',meta.definitions.windows{i}));
            end
            
            if i == 1
                ylabel(sprintf(freq_text{f}));
            end
            
            set(gca,'fontsize',20)
            
        end
    end
    filename = [name,'_spike'];
    print([plot_folder,filename],'-depsc');
    
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
               xlabel('Window') 
            end
            yticklabels([])
            
            if i == 1
                ylabel(sprintf(freq_text{f}));
            end
            set(gca,'fontsize',20)
        end
    end
    filename = [name,'_spike_network'];
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
        
    end
    %}
    yticklabels([])
    xlabel('Time (s)')
    title('Signal deviation')
    set(gca,'fontsize',20)
    
    axes(ha2(2))
    plot(avg_bin_dev,'ks-','linewidth',2)
    hold on
    title('Binned signal deviation')
    yticklabels([])
    xlabel('Window')
    set(gca,'fontsize',20)
    
    filename = [name,'_spike_dev'];
    print([plot_folder,filename],'-depsc');
    
end



end