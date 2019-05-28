function spike_plot(whichPts)

%% Parameters
sp_net = 1;

% 1 = alpha/theta; 2 = beta, 3 = low gamma, 4 = high gamma, 5 = ultra high, 6 = broadband
freq_text = {'alpha/theta','beta','low\ngamma','high\ngamma','ultra high\ngamma','broadband'};
%freq_text = {'alpha/theta'};
t_text = {'t1','t2','t3','t4','t5','t6','t7','t8','t9','t10','t11'};
n_f = length(freq_text);
n_times = 11*2+1;

%% Get file locations, load spike times and pt structure
if sp_net == 1
    addend_text = '_sp';
else
    addend_text = '';
end

locations = spike_network_files;
main_folder = locations.main_folder;
results_folder = [main_folder,'results/'];
data_folder = [main_folder,'data/'];
script_folder = locations.script_folder;
addpath(genpath(script_folder));
spike_times_file = [data_folder,'spike_times/times.mat'];
pt_file = [data_folder,'spike_structures/pt.mat'];
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
    stats_folder = [pt_folder,'stats/'];
    
    % Load the stats file
    out = load([stats_folder,'stats.mat']);
    stats = out.out;
    
    % Get variables
    dev = stats.signal.dev;
    bin_dev = stats.signal.bin_dev;
    ec = stats.network.ec;
    if sp_net == 0
        ge = stats.network.ge;
        ns = stats.network.ns;
        sync = stats.network.sync;
    elseif sp_net == 1
        ge = stats.network.ge_sp;
        ns = stats.network.ns_sp;
        sync = stats.network.sync;
    end
    
    
    avg_bin_dev = nanmean(bin_dev,1);

    plot_thing(1,:,:) = nanmean(ns,2);
    plot_thing(2,:,:) = nanmean(ec,2);
    plot_thing(3,:,:) = nanmean(ge,2);
    plot_thing(4,:,:) = nanmean(sync,2);

    orig_thing(1,:,:,:) = ns;
    orig_thing(2,:,:,:) = ec;
    orig_thing(3,:,:,:) = ge;
    orig_thing(4,:,:,:) = sync;

    plot_title{1} = 'Node strength\nof spike sequence chs';
    plot_title{2} = 'Eigenvector centrality\nof spike sequence chs';
    plot_title{3} = 'Global efficiency';
    plot_title{4} = 'Synchronizability';
    
    
    
    
    
    %% Compare ec spikes to ec not spikes (sanity check)
    %{
    if isfield(stats.network,'ec_notseq') == 1
        p = ranksum(squeeze(stats.network.ec_notseq(1,:,12)),...
            squeeze(stats.network.ec(1,:,12)));
        figure
        boxplot([squeeze(stats.network.ec_notseq(1,:,12))',squeeze(stats.network.ec(1,:,12))']);
        hold on
        plot([1 2],[0.2 0.2],'k-','linewidth',2)
        if p < 0.001
            text(1.5,0.2+0.01,sprintf('p < 0.001'),'fontsize',20,'horizontalalignment','center');
        else
            text(1.5,0.2+0.01,sprintf('p = %1.3f',p),'fontsize',20,'horizontalalignment','center');
        end
        xticklabels({'Non-spike electrodes','Spike electrodes'})
        ylabel('Mean eigenvector centrality');
        title('Eigenvector centrality during spike')
        set(gca,'fontsize',20);
        pause
        print([plot_folder,'ec_check'],'-depsc');
        close(gcf)
    end
    %}
    
    %% Plot aggregated metrics
    figure
    set(gcf,'position',[26 0 1242 900])
    [ha, pos] = tight_subplot(n_f-2, 4, [0.04 0.04], [0.08 0.08], [0.05 0.01]);
    for f = 1:n_f-2
        for i = 1:size(plot_thing,1)
            axes(ha((f-1)*4+i))
            plot(squeeze(plot_thing(i,f,:)),'ks-','linewidth',2)
            hold on
            for j = 2:size(plot_thing,3)
                h = ttest(squeeze(orig_thing(i,f,:,1)),...
                    squeeze(orig_thing(i,f,:,j)),'alpha',0.05/(length(avg_bin_dev)-1)/(n_f-2)/4);
                if h == 1
                    scatter(j,squeeze(plot_thing(i,f,j)),100,'r','filled')
                end
            end
            
            if f == 1
                title(sprintf(plot_title{i}));
            end
            if f == n_f-2
               xlabel('Time') 
            end
           % yticklabels([])
            
            if i == 1
                ylabel(sprintf(freq_text{f}));
            end
            set(gca,'fontsize',20)
        end
    end
    filename = [name,'_network_dev',addend_text];
    print([plot_folder,filename],'-depsc');
    
    figure
    set(gcf,'position',[26 0 600 700])
    [ha2, pos] = tight_subplot(2,1, [0.07 0.04], [0.09 0.06], [0.09 0.01]);
    axes(ha2(1))
    plot(avg_bin_dev,'ks-','linewidth',2)
    hold on
    for j = 2:length(avg_bin_dev)
        h = ttest(bin_dev(:,1),bin_dev(:,j),'alpha',0.05/(length(avg_bin_dev)-1));
        if h == 1
            scatter(j,avg_bin_dev(j),100,'r','filled');
        end
    end
    title('Binned signal z-score')
   % yticklabels([])
    %xlabel('Time (s)')
    xticklabels([])
    set(gca,'fontsize',20)
    
    axes(ha2(2))
    plot(avg_bin_dev,'ks-','linewidth',2)
    hold on
    for j = 2:length(avg_bin_dev)
        h = ttest(bin_dev(:,1),bin_dev(:,j),'alpha',0.05/(length(avg_bin_dev)-1));
        if h == 1
            scatter(j,avg_bin_dev(j),100,'r','filled');
        end
    end
    title('Binned signal deviation')
   % yticklabels([])
    xlabel('Time')
    set(gca,'fontsize',20)
    
    filename = [name,'_signal_dev_2',addend_text];
    print([plot_folder,filename],'-depsc');
   
end

end