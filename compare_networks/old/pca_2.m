function [ns,locs] = pca_2(whichPts)

which_freq = 2;

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
    
    % Skip it if name is empty
    if isempty(times(whichPt).name) == 1, continue; end
    name = times(whichPt).name;
    fprintf('\nDoing %s\n',name);
    pt_folder = [results_folder,name,'/'];
    adj_folder = [results_folder,name,'/just_spike/'];
    fs = pt(whichPt).fs;
    locs = pt(whichPt).new_elecs.locs;
    sz_times = pt(whichPt).newSzTimes;
    
    % Load thing
    stats_folder = [pt_folder,'stats/'];
    sp_adj = load([stats_folder,'sp_adj.mat']);
    sp_adj = sp_adj.sp_adj_info;
    
    adj = squeeze(sp_adj.adj(which_freq,:,:));
    
    % visualize it
    if 0
    figure
    imagesc(adj')
    end
    
    % Get avg adjacency matrix
    avg_adj = nanmean(squeeze(sp_adj.adj(which_freq,:,:)),1);
    avg_adj = flatten_or_expand_adj(avg_adj);
    ns = sum(avg_adj,1);
    
    if 0
        % plot the adjacency matrix
        figure
        subplot(1,2,1)
        imagesc(avg_adj)
        
        % Plot the node strength of each electrode
        subplot(1,2,2)
        scatter3(locs(:,1),locs(:,2),locs(:,3),100,'k')
        hold on
        scatter3(locs(:,1),locs(:,2),locs(:,3),100,sum(avg_adj,1),'filled')
        
    end
    
    
    
    % pca
    [coeff,score,latent] = pca(adj);
    
    % Plot the latencies
    if 0
    figure
    stem(latent)
    end
    
    % Get the first 2 subgraphs
    for i = 1:2
        
        % the flattened adj matrix for the principal component
        subgraph = coeff(:,i);
        
        % expand to get full adjacency matrix
        A = flatten_or_expand_adj(subgraph);
        nchs = size(A,1);
      
        
        
        if 0
        % Plot the histogram of scores of the component for the 1000
        % spikes
        subplot(1,2,2)
        histogram(score(:,i))
        end
        
        
        if 1
        % Show the adj and score for each electrode
        
        adj = coeff(:,i);
        adj = flatten_or_expand_adj(adj);
        
        figure
        set(gcf,'position',[135 364 1400 434])
        subplot(1,3,1)
        imagesc(adj)
        
        subplot(1,3,2)
        plot(sp_adj.time/3600,score(:,i))
        hold on
        for j = 1:size(sz_times,1)
            plot([sz_times(j,1) sz_times(j,1)]/3600,get(gca,'ylim'),'k')
        end
        
        
        subplot(1,3,3)
        scatter3(locs(:,1),locs(:,2),locs(:,3),100,'k');
        hold on
        scatter3(locs(:,1),locs(:,2),locs(:,3),100,sum(adj,1),'filled')
        end
        
        
        
    end
    
    
    
end


end