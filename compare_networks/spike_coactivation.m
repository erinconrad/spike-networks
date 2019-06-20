function spike_coactivation(whichPts)

%% File path
locations = spike_network_files;
spike_struct_folder = locations.spike_struct_folder;

%% Load files
pt = load([spike_struct_folder,'long_seq']);
pt = pt.pt;

cluster = load([spike_struct_folder,'cluster']);
cluster = cluster.cluster;

for whichPt = whichPts
    
    fprintf('Doing %s\n',pt(whichPt).name);
    
    %% Get patient parameters
    locs = pt(whichPt).electrodeData.locs(:,2:4);
    szTimes = pt(whichPt).newSzTimes;
    
    % Reorder seizure times if out of order
    oldSzTimes = szTimes;
    szTimes = sort(szTimes,1);
    if isequal(oldSzTimes,szTimes) == 0
        fprintf('WARNING!!! %s seizure times out of order\n',pt(whichPt).name);
    end
    
    % Combine nearly equal seizure times
    newIdx = 2;
    newSzTimes = [];
    newSzTimes(1,:) = szTimes(1,:);
    for j = 2:size(szTimes,1)
        if abs(szTimes(j,1)-szTimes(j-1,1)) < 10 && ...
                abs(szTimes(j,2)-szTimes(j-1,2))
           newIdx = newIdx - 1; 
           newSzTimes(newIdx,1) = min(szTimes(j,1),szTimes(j-1,1));
           newSzTimes(newIdx,2) = max(szTimes(j,2),szTimes(j-1,2));  
        else
           newSzTimes(newIdx,:) = szTimes(j,:);
        end
        newIdx = newIdx + 1;
    end
    
    if isequal(newSzTimes,szTimes) == 0
        fprintf('WARNING!!! %s had duplicate seizure times\n',pt(whichPt).name);
    end
    
    %% Get cluster info
    
    all_times_all = cluster(whichPt).all_times_all; % all spike times
    all_spikes = cluster(whichPt).all_spikes; % all spike channels
    all_locs = cluster(whichPt).all_locs;
    k = cluster(whichPt).k; % the number of clusters
    idx = cluster(whichPt).idx; % the cluster index for every spike
    C = cluster(whichPt).C; % the centroids of the clusters
    bad_cluster = cluster(whichPt).bad_cluster; % which clusters are bad
    
    % Confirm that I do not have any ictal spikes
    t = find(any(all_times_all >= szTimes(:,1)' & all_times_all <= szTimes(:,2)',2));
    if isempty(t) == 0
        fprintf('WARNING: Remaining ictal spikes for %s!\n',pt(whichPt).name);
        all_times_all(t) = [];
        all_spikes(t) = [];
        all_locs(t,:) = [];
        idx(t) = [];
    end
    

    % Remove bad clusters
    bad_idx = find(ismember(idx,bad_cluster));
    all_times_all(bad_idx) = [];
    all_spikes(bad_idx) = [];
    all_locs(bad_idx,:) = [];
    idx(bad_idx) = [];
    clusters = 1:k; clusters(bad_cluster) = [];
    C(bad_cluster,:) = [];
    popular = mode(idx); %most popular cluster
    nchs = size(locs,1);
    
    
    %% Initialize coA array
    window = 1800;
    time_thresh = 0.15;
    nbins = ceil((all_times_all(end)-all_times_all(1))/window);
    coA = zeros(nbins,nchs*(nchs-1)/2);
    times = zeros(nbins,1);
    counts = zeros(nbins,1);
    
    % Loop through and find coactivated channels
    for tt = 1:nbins
        
        % Get appropriate times
        curr_times = [all_times_all(1) + (tt-1)*window, min(all_times_all(1) + tt*window,all_times_all(end))];
        times(tt) = curr_times(2);
        
        % Get appropriate spikes
        sp_idx = find(all_times_all >= curr_times(1) & all_times_all <= curr_times(2));
        sp_times = all_times_all(sp_idx);
        sp_chs = all_spikes(sp_idx);
        
        coA_temp = zeros(nchs,nchs);
        
        counts(tt) = length(sp_times);
        
        % Loop through spikes
        for i = 1:length(sp_times)
            for j = i+1:length(sp_times)
                
                % if the two spikes are close enough, add to coactivation
                if sp_times(j) - sp_times(i) < time_thresh
                    coA_temp(sp_chs(i),sp_chs(j)) = coA_temp(sp_chs(i),sp_chs(j)) + 1;
                    coA_temp(sp_chs(j),sp_chs(i)) = coA_temp(sp_chs(j),sp_chs(i)) + 1;
                else
                    % if j is too far from i then no subsequent j will be
                    % closer, so need to move to the next i
                    break
                end
            end
        end
        
        % subtract mean
        coA_temp = coA_temp - mean(mean(coA_temp));
        
        % Flatten the adjacency matrix
        adj_flat = zeros(nchs*(nchs-1)/2,1);
        flat_idx = 0;

        for i = 1:nchs
            for j = 1:i-1
                flat_idx = flat_idx + 1;
                adj_flat(flat_idx) = coA_temp(j,i);
            end
        end
        
        
        % Add to output
        coA(tt,:) = adj_flat;
    end
    
    % Get avg adjacency matrix
    avg_adj = nanmean(coA,1);
    avg_adj = flatten_or_expand_adj(avg_adj);
    
    if 1 
       figure
       subplot(1,2,1)
       imagesc(avg_adj)
       
       subplot(1,2,2)
       scatter3(locs(:,1),locs(:,2),locs(:,3),100,'k')
       hold on
       scatter3(locs(:,1),locs(:,2),locs(:,3),100,sum(avg_adj,1),'filled')
    end
    
    
    % pca
    [coeff,score,latent] = pca(coA);
    
    if 0
        figure
        stem(latent)
    end
    
    if 0
    figure
    imagesc(coA');
    end
    
    if 0
    figure
    plot(times,counts)
    hold on
    for j = 1:size(szTimes,1)
        plot([szTimes(j,1) szTimes(j,1)],get(gca,'ylim'),'k')
    end
    end
    
    if 0
    % take first two principal components
    for i = 1:2
        adj = coeff(:,i);
        adj = flatten_or_expand_adj(adj);
        
        figure
        set(gcf,'position',[135 364 1400 434])
        subplot(1,3,1)
        imagesc(adj)
        subplot(1,3,2)
        plot(times,score(:,i))
        hold on
        for j = 1:size(szTimes,1)
            plot([szTimes(j,1) szTimes(j,1)],get(gca,'ylim'),'k')
        end
        subplot(1,3,3)
        scatter3(locs(:,1),locs(:,2),locs(:,3),100,'k')
        hold on
        scatter3(locs(:,1),locs(:,2),locs(:,3),100,sum(adj,1),'filled')

    end
    end
    
   
    
    
end


end