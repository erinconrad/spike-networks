function pca_2(whichPts)

which_freq = 2; % beta

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
    
    % Get number of channels
    y = size(adj,2);
    nchs = 0.5 + sqrt(0.25+2*y);
    
    if nchs-floor(nchs) > 0.01
        error('what\n');
    end
    
    nchs = round(nchs);
    
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
        
        % Reconstruct the full adjancency matrix
        A = zeros(nchs,nchs);
        count = 0;
        for j = 1:size(A,1)
            for k = 1:j-1
                count = count + 1;
                A(j,k) = subgraph(count);
            end
        end
        
        if count ~= length(subgraph)
            error('what\n');
        end
        
        A = A + A';
        
        % Plot the adjacency matrix for the first principal component
        figure
        subplot(1,2,1)
        imagesc(A)
        
        % Plot the histogram of scores of the first component for the 1000
        % spikes
        subplot(1,2,2)
        histogram(score(:,i))
        
        % Find the 50 spikes with the highest score
        [~,sorted_score_I] = sort(score(:,i));
        top_score = sorted_score_I(end-49:end);
        bottom_score = sorted_score_I(1:50);
        
        fprintf('Top score for component %d:\n',i);
        sp_adj.labels(top_score)
        
        fprintf('Bottom score for component %d:\n',i);
        sp_adj.labels(bottom_score)
        
        fprintf('\n');
        
    end
    
    
end


end