function compare_spikes_by_loc(main,spike,pre,timing,pts)

%% Locations
locations = spike_network_files;
main_folder = locations.main_folder;
data_folder = [main_folder,'data/'];
pt_folder = [data_folder,'spike_structures/'];

%% Load pt structure
pt = load([pt_folder,'pt.mat']);
pt = pt.pt;

npts = length(pts);
all_change = cell(npts,1);
all_rel_change = cell(npts,1);

for i = 1:npts
    p = pts(i);

    %% Get max spike power electrode
    sp_power = spike(p).spike_powers;
    [~,max_ch] = max(sp_power,[],2);
    
    %% Get electrode locs
    locs = [];
    name = main.pt(p).name;
    for j = 1:length(pt)
        if strcmp(pt(j).name,name)
            locs = pt(j).new_elecs.locs;
        end
    end
    
    %% initialize arrays
    % Get sizes of things
    sp = 1;
    data = main.pt(p).sp_or_not(sp).data;
    main_ch_data = nan(size(data,1),size(data,2));
    first = nan(size(data,1),1);
    last = nan(size(data,1),1);

    for s = 1:size(data,1)
        % Get the power in the biggest spike power electrode
        main_ch_data(s,:) = data(s,:,max_ch(s));
        current_pre = pre(p).before_rise(s,:);
        % remove non-pre times
        last_pre = find(current_pre ==0);
        last_pre = last_pre(1) - 1; % the last pre-spike window

        first(s) = main_ch_data(s,1);
        last(s) = main_ch_data(s,last_pre);
    end

    %% Get change
    change = last-first;
    all_change{i} = change;
    all_rel_change{i} = (last-first)./abs(first);
    
    %
    %% Only include most common spikes
    % Get the n channels corresponding to 80% of all spikes for that
    % patient (to ignore the one-off channels that will reduce my power)
    [C,ia,ic] = unique(max_ch); % get unique values and occurrences
    prct_spike = 0.7;
    num_occurrences = [C, accumarray(ic, 1)];
    num_sorted = sortrows(num_occurrences,2,'descend');
    perc_captured = cumsum(num_sorted(:,2))/length(change);
    first_above_thresh = find(perc_captured > prct_spike);
    if isempty(first_above_thresh)
        first_above_thresh = size(num_sorted,1);
    else
        first_above_thresh = first_above_thresh(1);
    end
    chs_above_thresh = num_sorted(1:first_above_thresh,1);
    bin_sp = ismember(max_ch,chs_above_thresh);
    %}

    %% Cluster by location and compare change between clusters
    %{
    eva = evalclusters(locs(max_ch,:),'kmeans','CalinskiHarabasz','KList',1:length(unique(max_ch)));
    k = eva.OptimalK;
    idx = kmeans(locs(max_ch,:),3);
    %}
    
    %
    sumds = nan(length(unique(max_ch)),1);
    for k = 1:length(unique(max_ch))
        [~,~,sumd] = kmeans(locs(max_ch,:),k);
        sumds(k) = sum(sumd);
    end
    if 1
        figure
        plot(sumds);
        fprintf('\nSelect the cluster number\n');
        [x,~] = ginput;
        close(gcf)
    end
    k = round(x);
    idx = kmeans(locs(max_ch,:),k);
    %}
    
    if 1
    figure
    set(gcf,'position',[71 400 1370 398]);
    subplot(1,2,1)
    scatter3(locs(:,1),locs(:,2),locs(:,3),100,'o')
    hold on
    scatter3(locs(max_ch,1),locs(max_ch,2),locs(max_ch,3),100,idx,'o','filled')
    
    % anova, where the group is the cluster index
    [pval,tbl] = anova1(change,idx,'off');
    subplot(1,2,2)
    scatter(idx,change)
    title(sprintf('anova p-value: %1.3f',pval))
    pause
    close(gcf)
    end
    
    %% Compare change between clusters
    
    
    %
    % Categorical way, ignoring location
    
    
    if 0
    figure
    [pval,tbl] = anova1(change(bin_sp),max_ch(bin_sp),'off');
    plot(max_ch(bin_sp),change(bin_sp),'o')
    hold on
    plot(xlim,[0 0],'k--');
    title(sprintf('Unique channels: %d\nOne way anova p value %1.3f',...
        length(unique(max_ch(bin_sp))),pval));
    % one way anova
    
    pause
    close(gcf)
    end
    %}
    
end


end