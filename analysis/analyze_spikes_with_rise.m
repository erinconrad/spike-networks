function all_rel_change = analyze_spikes_with_rise(main,spike,pre,timing,pts)

change_thresh = 5;
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
    
    
    %% Does spike location determine pre-spike rise?
    %{
    Seems like no, more or less randomly distributed
    %}
    if 0
    [pval,tbl] = anova1(change,max_ch,'off');
    plot(max_ch,change,'o')
    title(sprintf('One way anova p value %1.3f',pval));
    % one way anova
    
    pause
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