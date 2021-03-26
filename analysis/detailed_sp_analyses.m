function detailed_sp_analyses(rel_change,rise_pts,spike_power,pre_spike,abs_power,ers,ns)

alpha = 0.05;
nfreq = 3;

all_earliest = nan(length(rise_pts),1);

%% When does the change occur?
% Loop over rise pts
for i = 1:length(rise_pts)
    
    p = rise_pts(i);

    %% Get max spike power electrode
    sp_power = spike_power(p).spike_powers;
    [~,max_ch] = max(sp_power,[],2);

    %% Get before rise data
    before_rise = pre_spike(p).before_rise;

    %% Find those spikes with a large relative power change
    %{
    curr_pt_change = rel_change{i};
    change_sp_idx = find(curr_pt_change > rel_change_thresh);
    %}
    nsp = size(abs_power.freq(1).pt(p).sp_or_not(1).data,1);
    
    all_data_before_rise = nan(nsp,sum(before_rise(1,:)));
    
    % Loop over spikes with change
    for s = 1:nsp%length(change_sp_idx)
        
        % get spike and max ch
        %s = change_sp_idx(j);
        curr_max_ch = max_ch(s);
        
        data = abs_power.freq(1).pt(p).sp_or_not(1).data;

        % spike data for max channel over time
        sp_data = squeeze(data(s,:,curr_max_ch));
        sp_before_rise = logical(before_rise(s,:));
        
        % just before rise times
        sp_data_before_rise = sp_data(sp_before_rise);
        all_data_before_rise(s,1:length(sp_data_before_rise)) = sp_data_before_rise;

    end
    
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
    for t = 2:size(all_data_before_rise,2)
        num_left = size(all_data_before_rise,2) - t + 1; % how many left including this one
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

all_earliest
        
%% What frequencies are involved in the change
for f = 1:nfreq
    curr_freq = ers.freq(f);
    agg_array = nan(length(rise_pts),1);
    pval_array = nan(length(rise_pts),1);
    % Loop over pts
    for i = 1:length(rise_pts)
        p = rise_pts(i);
        data = curr_freq.pt(p).sp_or_not(1).data;
        nsp = size(data,1);
        
        % Get max spike power electrode
        sp_power = spike_power(p).spike_powers;
        [~,max_ch] = max(sp_power,[],2);
        
        % pre-rise data
        before_rise = pre_spike(p).before_rise;
        all_change = nan(nsp,2);
        
        for s = 1:nsp
            
            curr_max_ch = max_ch(s);
            sp_data = squeeze(data(s,:,curr_max_ch));
            first = sp_data(1);
            
            current_pre = before_rise(s,:);
            last_pre = find(current_pre==0);
            last_pre = last_pre(1) - 1;
            last = sp_data(last_pre);
            
            all_change(s,:) = [first,last];
            
        end
        
        if 0
            plot(1+0.05*rand(size(all_change,1),1),all_change(:,1),'o')
            hold on
            plot(2+0.05*rand(size(all_change,1),1),all_change(:,2),'o')
        end
        
        [~,pval,~,stats] = ttest(all_change(:,1),all_change(:,2));
        pval_array(i) = pval;
        agg_array(i) = stats.tstat;
    end
    
    %pval_array
    [~,pval] = ttest(agg_array);
    fprintf('For freq %d, aggregate pval is %1.3f\n',f,pval);
end

%% What electrodes are involved

%% Is there a corresponding network change?
    







end


