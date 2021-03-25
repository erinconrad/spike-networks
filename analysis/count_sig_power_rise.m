function count_sig_power_rise(main,spike,clean_pre,f,alpha)

np = length(spike);

% Initialize p-value array
pval_array = nan(np,2);
agg_array = nan(np,2);

% Loop over patients 
for p = 1:np
    
    %% Get max spike power electrode
    sp_power = spike(p).spike_powers;
    [~,max_ch] = max(sp_power,[],2);

    %% Get fake max spike power electrode for non-spike periods
    n_not_spikes = size(main.freq(f).pt(p).sp_or_not(2).data,1);
    % Randomly sample n_not_spikes from max chs
    fake_max_ch = datasample(max_ch,n_not_spikes);

    %% Get fake pre-times for non-spike periods
    for s = 1:n_not_spikes
        clean_pre(p).fake_before_rise(s,:) = mode(clean_pre(p).before_rise,1);
    end
    
    for sp = 1:2
        data = main.freq(f).pt(p).sp_or_not(sp).data;
        main_ch_data = nan(size(data,1),size(data,2));
        first = nan(size(data,1),1);
        last = nan(size(data,1),1);

        % Loop over spikes
        for s = 1:size(data,1)

            if sp == 1
                % Get the power in the biggest spike power electrode
                main_ch_data(s,:) = data(s,:,max_ch(s));
                current_pre = clean_pre(p).before_rise(s,:);
            else
                % Get the power in the biggest fake power electrode
                main_ch_data(s,:) = data(s,:,fake_max_ch(s));
                current_pre = clean_pre(p).fake_before_rise(s,:);
            end

            % remove non-pre times
            last_pre = find(current_pre ==0);
            last_pre = last_pre(1) - 1; % the last pre-spike window
            main_ch_data(s,current_pre==0) = nan;
            

            first(s) = main_ch_data(s,1);
            last(s) = main_ch_data(s,last_pre);

        end
        
        % Paired t-test comparing first and last time
        [~,pval,~,stats] = ttest(last,first);
        pval_array(p,sp) = pval;
        
        agg_array(p,sp) = stats.tstat;
        
    end
    
end

fprintf('\n%d of %d patients had a significant pre-IED rise.\n',...
    sum(pval_array(:,1)<alpha),size(pval_array,1));

fprintf('\n%d of %d patients had a significant rise in IED-free periods.\n',...
    sum(pval_array(:,2)<alpha),size(pval_array,1));

% Do aggregate test
[~,pval_spike] = ttest(agg_array(:,1));
[~,pval_not] = ttest(agg_array(:,2));

fprintf('\nThe aggregate p-vals are %1.3f for IED and %1.3f for IED-free.\n',...
    pval_spike,pval_not);

end