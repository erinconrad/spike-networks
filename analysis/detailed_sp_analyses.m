function detailed_sp_analyses(rise_pts,spike_power,pre_spike,abs_power,ers,ns)


alpha = 0.05;
nfreq = 3;


%% When does the change occur?
rel_change_thresh = -inf;
all_earliest = nan(length(rise_pts),1);
% Loop over rise pts
for i = 1:length(rise_pts)
    
    p = rise_pts(i);

    %% Get max spike power electrode
    sp_power = spike_power(p).spike_powers;
    [~,max_ch] = max(sp_power,[],2);

    %% Get before rise data
    before_rise = pre_spike(p).before_rise;
    
    % Get last time before any rise
    any_before_rise = sum(any(before_rise,1));

    %% Find those spikes with a large relative power change
    % Get spike abs power data
    data = abs_power.freq(1).pt(p).sp_or_not(1).data;
    
    % get the data on the channel with the highest spike power change
    data_main_ch = get_data_from_spec_ch(data,max_ch);
    
    % calcaulte the relative power change compared to baseline
    rel_change = calc_rel_change(data_main_ch,before_rise,1:5);
    
    % decide if it's above a threshold
    change_sp_idx = find(rel_change > rel_change_thresh);
    
    %nsp = size(abs_power.freq(1).pt(p).sp_or_not(1).data,1);
    

    % get the data occurring before the rise (setting all other times to
    % nans)
    data_change_sp_idx = data_main_ch(change_sp_idx,:);
    before_rise_change_sp_idx = before_rise(change_sp_idx,:);
    all_data_before_rise = data_change_sp_idx;
    all_data_before_rise(before_rise_change_sp_idx==0) = nan; % set times after the rise to nans
    
    
    time_pvals = nan(size(all_data_before_rise,2),1);
    % Loop over times to see which have significant difference from first
    % time
    for t = 2:size(all_data_before_rise,2)
        
        % Do paired ttest
        [~,pval] = ttest(all_data_before_rise(:,1),all_data_before_rise(:,t));
        time_pvals(t) = pval;
    end

    % See which have p < 0.05 for it and all subsequent ones
    earliest = nan;
    less_alpha = time_pvals < alpha;
    for t = 2:any_before_rise
        num_left = any_before_rise - t + 1; % how many left including this one
        num_sub_alpha = sum(less_alpha(t:length(less_alpha)));
        if num_left == num_sub_alpha
            earliest = t;
            break
        end
    end
    
    if isnan(earliest)
        earliest_time = nan;
    else
        earliest_time = abs_power.times(earliest);
    end
    
    all_earliest(i) = earliest_time;
end

% Display earliest time of rise
all_earliest
        
%% What frequencies are involved in the change

change_by_f = cell(length(rise_pts),1);
for i = 1:length(rise_pts)
    p = rise_pts(i);
    change_by_f{i} = nan(size(abs_power.freq(1).pt(p).sp_or_not(1).data,1),nfreq);
end
for f = 1:nfreq
    for i = 1:length(rise_pts)
        
        p = rise_pts(i);

        %% Get max spike power electrode
        sp_power = spike_power(p).spike_powers;
        [~,max_ch] = max(sp_power,[],2);

        %% Get before rise data
        before_rise = pre_spike(p).before_rise;

        %% Find those spikes with a large relative power change
        % Get spike abs power data
        data = abs_power.freq(1).pt(p).sp_or_not(1).data;

        % get the data on the channel with the highest spike power change
        data_main_ch = get_data_from_spec_ch(data,max_ch);

        % calcaulte the relative power change compared to baseline
        rel_change = calc_rel_change(data_main_ch,before_rise,1:5);

        % decide if it's above a threshold
        change_sp_idx = find(rel_change > rel_change_thresh);

        % get the data occurring before the rise (setting all other times to
        % nans)
        data = ers.freq(f).pt(p).sp_or_not(1).data;
        data_main_ch = get_data_from_spec_ch(data,max_ch);
        data_change_sp_idx = data_main_ch(change_sp_idx,:);
        before_rise_change_sp_idx = before_rise(change_sp_idx,:);
        all_data_before_rise = data_change_sp_idx;
        all_data_before_rise(before_rise_change_sp_idx==0) = nan; % set times after the rise to nans
    
        % loop through each spike and get first and last
        for s = 1:size(all_data_before_rise,1)
            curr = all_data_before_rise(s,:);
            first = curr(1);
            last_idx = find(isnan(curr));
            last_idx = last_idx(1) - 1;
            last = curr(last_idx);
            change_by_f{i}(s,f) = last-first;
        end
        
        
        
        
    end

           
end


% compare with anova - within patient
mean_changes = nan(length(rise_pts),nfreq);
for i = 1:length(rise_pts)
    changes = change_by_f{i};
    pval = anova1(changes,[],'off');
    
    if 0
        plot(1+0.05*rand(size(change_by_f{i},1),1),change_by_f{i}(:,1),'o')
        hold on
        plot(2+0.05*rand(size(change_by_f{i},1),1),change_by_f{i}(:,2),'o')
        plot(3+0.05*rand(size(change_by_f{i},1),1),change_by_f{i}(:,3),'o')
        title(sprintf('sub-gamma: %1.1e, low-gamma: %1.1e, high-gamma: %1.1e\np = %1.3f',...
            mean(changes(:,1)),mean(changes(:,2)),mean(changes(:,3)),pval))
        pause
        hold off
    end
    
    for f = 1:nfreq
        mean_changes(i,f) = mean(changes(:,f));
    end
    
end

% Compare frequency changes across patients
pval = anova1(mean_changes,[],'off');
if 0
    for f = 1:nfreq
        plot(f+0.05*rand(size(mean_changes,1),1),mean_changes(:,f),'o')
        hold on
        title(sprintf('sub-gamma: %1.1e, low-gamma: %1.1e, high-gamma: %1.1e\np = %1.3f',...
            mean(mean_changes(:,1)),mean(mean_changes(:,2)),mean(mean_changes(:,3)),pval))
    end
end

%% What electrodes are involved
all_rho_and_p = nan(length(rise_pts),2);
% Loop over rise pts
for i = 1:length(rise_pts)
    p = rise_pts(i);

    
    % Get power in spike periods across channels
    sp_power = spike_power(p).spike_powers;
    
    % Get pre-spike abs power data
    data = abs_power.freq(1).pt(p).sp_or_not(1).data;
    
    % before rise indicators
    before_rise = pre_spike(p).before_rise;

    % Loop over spikes
    nsp = size(data,1);
    nch = size(data,3);
    sp_rhos = nan(nsp,1);
    sp_change = nan(nsp,nch);
    sp_rel_change = nan(nsp,nch);
    for s = 1:nsp
        curr_sp_power = squeeze(sp_power(s,:));
        curr_sp_pre_data = squeeze(data(s,:,:));
        curr_before_rise = squeeze(before_rise(s,:));
        
        last_pre_rise = find(curr_before_rise == 0);
        last_pre_rise = last_pre_rise(1) - 1; % move one back to get last pre rise
        
        curr_sp_first = curr_sp_pre_data(1,:);
        curr_sp_last = curr_sp_pre_data(last_pre_rise,:);
        change = (curr_sp_last - curr_sp_first);
        rel_change = (curr_sp_last - curr_sp_first)./abs(curr_sp_first);
        
        curr_rho = corr(curr_sp_power',change');
        sp_rhos(s) = curr_rho;
        sp_change(s,:) = change;
        sp_rel_change(s,:) = rel_change;
    end
    
    % Average change
    mean_sp_change = mean(sp_change,1);
    mean_sp_rel_change = mean(sp_rel_change,1);
    mean_sp_power = mean(sp_power,1);
    
    [rho,pval] = corr(mean_sp_power',mean_sp_change');
    all_rho_and_p(i,:) = [rho,pval];
    
end



%% Is there a corresponding network change?
% Loop over rise pts
for i = 1:length(rise_pts)
    p = rise_pts(i);
    
    % Get power in spike periods across channels
    sp_power = spike_power(p).spike_powers;
    [~,max_ch] = max(sp_power,[],2);
    
    % before rise indicators
    before_rise = pre_spike(p).before_rise;
    nsp = size(sp_power,1);
    
    first_last = nan(nsp,nfreq,2);
    
    for f = 1:nfreq
    
        % Get node strength
        data = ns.freq(f).pt(p).sp_or_not(1).data;

        % get the data on the channel with the highest spike power change
        data_main_ch = get_data_from_spec_ch(data,max_ch);
    
        nsp = size(data_main_ch,1);
        
        
        % Loop over spikes
        for s = 1:nsp
            data_sp = squeeze(data_main_ch(s,:));
            before_rise_sp = squeeze(before_rise(s,:));
            last_time = find(before_rise_sp==0);
            last_time = last_time(1)-1; % 1 before is the last before rise
            
            first = data_sp(1);
            last = data_sp(last_time);
            first_last(s,f,:) = [first last];
        end
        
    end
    
    % Plot and significance testing on single pt level
    if 1
        figure
        set(gcf,'position',[440 490 1001 308])
        for f = 1:nfreq
            [~,pval] = ttest(squeeze(first_last(:,f,2)),squeeze(first_last(:,f,1)));
            
            plot(f+0.05*rand(size(first_last,1),1),...
                first_last(:,f,2)-first_last(:,f,1),'o') 
            hold on
            text(f,max(first_last(:,f,2)-first_last(:,f,1)),...
                sprintf('%1.3f',pval),'fontsize',20)
        end
        plot(xlim,[0 0],'k')
        pause
        close(gcf)
    end
    
end







end


