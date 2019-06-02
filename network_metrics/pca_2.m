function pca_2(whichPts)

which_freq = 4; % beta

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
        
        % Plot the adjacency matrix for the first principal component
        figure
        subplot(1,2,1)
        imagesc(A)
        
        
        % Plot the histogram of scores of the first component for the 1000
        % spikes
        subplot(1,2,2)
        histogram(score(:,i))
        
        % Find the 10 spikes with the highest score
        [~,sorted_score_I] = sort(score(:,i));
        top_score = sorted_score_I(end-9:end);
        bottom_score = sorted_score_I(1:10);
        
        fprintf('Top score for component %d:\n',i);
        sp_adj.labels(top_score)
        
        fprintf('Bottom score for component %d:\n',i);
        sp_adj.labels(bottom_score)
        
        fprintf('\n');
        
    end
    
    % Plot the scores over time
    figure
    for i = 1:3
        plot(sp_adj.time,score(:,i))
        hold on
    end
    
    
end


end