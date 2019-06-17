function pca_2(whichPts)

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
    
    
    
    % pca
    [coeff,score,latent] = pca(adj);
    
    % Plot the latencies
    if 0
    figure
    stem(latent)
    end
    
    % Get the first 3 subgraphs
    for i = 1:3
        
        % the flattened adj matrix for the first principal component
        subgraph = coeff(:,i);
        
        % expand to get full adjacency matrix
        A = flatten_or_expand_adj(subgraph);
        nchs = size(A,1);
        
        if 1
        % Plot the adjacency matrix for the first principal component
        figure
        subplot(1,2,1)
        imagesc(A)
        
        
        % Plot the histogram of scores of the first component for the 1000
        % spikes
        subplot(1,2,2)
        histogram(score(:,i))
        
        % Find the 50 spikes with the highest score and show the sequences
        % on a brain
        [~,sorted_score_I] = sort(score(:,i));
        top_score = sorted_score_I(end-49:end);
        bottom_score = sorted_score_I(1:50);
        
        locs = pt(whichPt).new_elecs.locs;
        figure
        subplot(1,2,1)
        scatter3(locs(:,1),locs(:,2),locs(:,3),100);
        hold on
        ch_count = zeros(nchs,1);
        for j = 1:length(top_score)
            if top_score(j) > length(sp_adj.labels), continue; end
            labels = sp_adj.labels(top_score(j)).labels;
            ch_nums = get_ch_nums_from_labels(pt,whichPt,labels);
            ch_nums(ch_nums==0) = [];
            ch_count(ch_nums) = ch_count(ch_nums) + 1;
        end
        
        scatter3(locs(:,1),locs(:,2),locs(:,3),100,ch_count,'filled');
        title(sprintf('Top scores for %d',i))
        
        subplot(1,2,2)
        
        scatter3(locs(:,1),locs(:,2),locs(:,3),100);
        hold on
        ch_count = zeros(nchs,1);
        for j = 1:length(bottom_score)
            labels = sp_adj.labels(bottom_score(j)).labels;
            ch_nums = get_ch_nums_from_labels(pt,whichPt,labels);
            ch_nums(ch_nums==0) = [];
            ch_count(ch_nums) = ch_count(ch_nums) + 1;
        end
        
        scatter3(locs(:,1),locs(:,2),locs(:,3),100,ch_count,'filled');
        title(sprintf('Bottom scores for %d',i))
        end
        
    end
    
    % Plot the scores over time
    figure
    for i = 1
        plot(sp_adj.time/3600,score(:,i))
        hold on
        for j = 1:size(sz_times,1)
            plot([sz_times(j,1) sz_times(j,1)]/3600,get(gca,'ylim'),'k')
        end
    end
    
    % Pretend the signal is just the first principal component
    new_adj_flat = coeff(:,1)*score(:,1)';
    
end


end