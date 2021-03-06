function calc_network_metrics(whichPts,small)

%% Parameters
% 1 = alpha/theta; 2 = beta, 3 = low gamma, 4 = high gamma, 5 = ultra high, 6 = broadband
freq_text = {'alpha/theta','beta','low\ngamma','high\ngamma','ultra high\ngamma','broadband'};
%freq_text = {'alpha/theta'};
n_f = length(freq_text);


%% Get file locations, load spike times and pt structure
locations = spike_network_files;
main_folder = locations.main_folder;
results_folder = [main_folder,'results/'];
data_folder = [main_folder,'data/'];
script_folder = locations.script_folder;
addpath(genpath(script_folder));
addpath(genpath(locations.BCT));
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
    if small == 1
        adj_folder = [results_folder,name,'/adj_small/'];
        stats_folder = [pt_folder,'stats_small/'];
    elseif small == 0
        adj_folder = [results_folder,name,'/adj/'];
        stats_folder = [pt_folder,'stats/'];
    elseif small == 2
        adj_folder = [results_folder,name,'/adj_test/'];
        stats_folder = [pt_folder,'stats_test/'];
    elseif small == 3
        adj_folder = [results_folder,name,'/adj_simple/'];
        stats_folder = [pt_folder,'stats_simple/'];
        n_f = 1;
    elseif small == 4
        adj_folder = [results_folder,name,'/adj_coherence/'];
        stats_folder = [pt_folder,'stats_coherence/'];
    elseif small == 5
        adj_folder = 'adj_mat/manual/adj_simple/';
        stats_folder = 'stats/manual/simple/';
    end
    fs = pt(whichPt).fs;
    %error('look\n');
    
    if exist(stats_folder,'dir') == 0
        mkdir(stats_folder);
    end
    
    %% Load adjacency matrices and calculate metrics
    
    if small == 5
        listing = dir([adj_folder,name,'*']);
    else
        listing = dir([adj_folder,'adj*.mat']);
        
    end
    
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
    
    
    n_times = size(meta.spike(1).index_windows,1);
    n_seconds = 14;
    
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
    ge_sp_net = nan(n_f,n_spikes,n_times);
    sync_sp_net = nan(n_f,n_spikes,n_times);
    ns_sp_net = nan(n_f,n_spikes,n_times);
    
    % Initialize spike count
    s_count = 0;
    out = [];
    
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
            index_windows = meta.spike(s).index_windows;
            
            % readjust size
            
            if 1
                if length(values) < size(dev,2)
                    fprintf('Padding values by %d\n',size(dev,2)-length(values));
                    values = [values;nan(size(dev,2)-length(values),1)];

                elseif length(values) > size(dev,2)
                    fprintf('Shortening values by %d\n',size(dev,2)-length(values));
                    values = values(1:end-(length(values)-size(dev,2)));

                end
            
                
            end
            
            dev_t = sqrt((values-nanmedian(values)).^2);
            dev(s_count,:) = dev_t;
            
            for i = 1:size(index_windows,1)
                bin_dev(s_count,i) = nanmean(dev_t(round(index_windows(i,1)):...
                    round(index_windows(i,2))));
            end
            
            for which_freq = 1:n_f
                if small == 3
                    adj_all_t= meta.spike(s).adj;
                else
                    adj_all_t= meta.spike(s).adj(which_freq).adj;
                    
                end

                % Loop through times
                for tt = 1:size(adj_all_t,1)

                    % Get adj matrix of interest
                    adj = squeeze(adj_all_t(tt,:,:));
                    
                    if sum(sum(isnan(adj))) > 0
                        continue
                    end
                    
                    % also get the adjacency matrix for the network made up
                    % of only the spike sequence channels
                    adj_sp = adj(seq_chs,seq_chs);
                    
                    if 0
                       % Side by side plot
                       figure
                       subplot(1,2,1)
                       imagesc(adj)
                       
                       subplot(1,2,2)
                       imagesc(adj_sp)
                       pause
                       close(gcf)
                       
                    end
                    
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
                    
                    %% Calculate metrics for network of just spike channels
                    if sum(seq_chs) > 4
                        ge_sp_net(which_freq,s_count,tt) = efficiency_wei(adj_sp,0);
                        sync_sp_net(which_freq,s_count,tt) = synchronizability_sp(adj_sp);

                        % mean node strength of all channels in spike network
                        ns_sp_net(which_freq,s_count,tt) = mean(strengths_und(adj_sp));
                    end

                end
                
            end
            
            out.index_windows(s_count).window = index_windows;
            
        end
        
    end
    
    %% Aggregate metrics
    %{
    avg_ns_seq = nanmean(ns_seq-median(ns_seq,3),2);
    avg_ge = nanmean(ge - median(ge,3),2);
    avg_sync = nanmean(sync - median(sync,3),2);
    avg_ec_seq = nanmean(ec_seq - median(ec_seq,3),2);
    %}
    
    %% Save stuff to structure
    
    
    out.fs = fs;
    out.signal.dev = dev;
    out.signal.bin_dev = bin_dev;
    
    out.network.ns = ns_seq;
    out.network.ec = ec_seq;
    out.network.ns_notseq = ns_not_seq;
    out.network.ec_notseq = ec_not_seq;
    out.network.ge = ge;
    out.network.sync = sync;
    out.network.ns_sp = ns_sp_net;
    out.network.ge_sp = ge_sp_net;
    out.network.sync_sp = sync_sp_net;
    out.network.index_windows = meta.spike(1).index_windows;

    if small == 1
        save([stats_folder,'stats_small.mat'],'out');
    elseif small == 0
        save([stats_folder,'stats.mat'],'out');
    elseif small == 2
        save([stats_folder,'stats_test.mat'],'out');
    elseif small == 3
        save([stats_folder,'stats_simple.mat'],'out');
    elseif small == 4
        save([stats_folder,'stats_coherence.mat'],'out');
    end
    
    z_ns = (((ns_seq-mean(ns_seq,3))./std(ns_seq,0,3)));
    z_ec = (((ec_seq-mean(ec_seq,3))./std(ec_seq,0,3)));
    z_sync = (((sync-mean(sync,3))./std(sync,0,3)));
    z_ge = (((ge-mean(ge,3))./std(ge,0,3)));
    
    avg_z_ns = nanmean(z_ns,2);
    avg_z_ec = nanmean(z_ec,2);
    avg_z_sync = nanmean(z_sync,2);
    avg_z_ge = nanmean(z_ge,2);
    
    plot_thing(1,:,:) = avg_z_ns;
    plot_thing(2,:,:) = avg_z_ec;
    plot_thing(3,:,:) = avg_z_ge;
    plot_thing(4,:,:) = avg_z_sync;
    
    plot_title{1} = 'Node strength\nof spike sequence chs';
    plot_title{2} = 'Eigenvector centrality\nof spike sequence chs';
    plot_title{3} = 'Global efficiency';
    plot_title{4} = 'Synchronizability';
    
    %% Get avg deviation of signal
    
    avg_dev = nanmean(dev,1);
    avg_bin_dev = nanmean(bin_dev,1);
    %}
    

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
               xlabel('Window') 
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
    
    filename = [name,'_signal_dev'];
    print([plot_folder,filename],'-depsc');
    
   
    
end


end