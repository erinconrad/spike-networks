function all_rel_change = analyze_spikes_with_rise(main,spike,pre,timing,pts)

change_thresh = 5;


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

    %% Get times
    sp_times = timing(p).times;
    
    % Get electrode locs
    locs = [];
    name = main.pt(p).name;
    for j = 1:length(pt)
        if strcmp(pt(j).name,name)
            locs = pt(j).new_elecs.locs;
        end
    end

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
    
    %% Cluster spikes by location
    
    % Get locations of the max ch for each spike
    spike_locs = locs(max_ch,:);
    
    
    %% Does spike electrode determine pre-spike rise? (Basic categorical test)
    prct_spike = 1;
    
    % Get the n channels corresponding to 80% of all spikes for that
    % patient (to ignore the one-off channels that will reduce my power)
    [C,ia,ic] = unique(max_ch); % get unique values and occurrences
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
    
    if 0
    figure
    [pval,tbl] = anova1(change(bin_sp),max_ch(bin_sp),'off');
    plot(max_ch(bin_sp),change(bin_sp),'o')
    hold on
    plot(xlim,[0 0],'k--');
    title(sprintf('One way anova p value %1.3f',pval));
    % one way anova
    
    pause
    close(gcf)
    end
    
    %% Does spike electrode location determine pre-spike rise? 
    if 0
        mean_change = nan(size(locs,1),1);
        for ich = 1:size(locs,1)
            spikes_on_ch_bin = max_ch == ich;
            change_on_ch = change(spikes_on_ch_bin);
            mean_change(ich) = mean(change_on_ch);
        end
        wij = getwij(locs,20); % adjust d
        MI = moranStats(mean_change',wij);
        
        figure

        for ich = 1:size(locs,1)
            
            scatter3(locs(ich,1),locs(ich,2),locs(ich,3),100,...
                mean_change(ich),'filled');
            hold on
        end
        scatter3(locs(:,1),locs(:,2),locs(:,3),100,'k');
        colorbar
        title(sprintf('Moran index: %1.2f, p = %1.3f',MI.I,MI.p));
        pause
        close(gcf)
    end
    
    %% Do spikes with pre-spike changes tend to occur at different times
    %{
    Seems like no, more or less randomly distributed
    %}
    if 0
        % dumb way
        [r,pval] = corr(sp_times,change);
        if isnan(r), error('what'); end
        plot(sp_times,change,'o')
        title(sprintf('Pearson corr r = %1.2f, p = %1.3f',r,pval))
        pause
    end

end








end