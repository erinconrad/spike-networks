function metrics = generate_summary_stats(metrics,met,include_times,rm_rise,...
    is_spike_soz,do_cumulative,is_spike_depth)

nfreq = length(metrics.time.freq);

for f = 1:nfreq
    
    
    % Prep array
    metrics.time.freq(f).(met).auc.data = zeros(...
        length(metrics.time.freq(f).(met).pt),2);  
    metrics.time.freq(f).(met).first_last = zeros(...
        length(metrics.time.freq(f).(met).pt),1);
    metrics.time.freq(f).(met).auc.individual_ranksum_p = zeros(...
        length(metrics.time.freq(f).(met).pt),1);
    metrics.time.freq(f).(met).auc.info = {'spike','not'};
    metrics.time.freq(f).(met).median_data.data = zeros(...
        length(metrics.time.freq(f).(met).pt),...
        size(metrics.time.freq(f).(met).pt(1).spike.data,2),2); 
    metrics.time.freq(f).(met).median_data.info = {'spike','not'};
    metrics.time.freq(f).(met).auc.times = metrics.time.freq(f).(met).pt(1).times;
    metrics.time.freq(f).(met).auc.short.data = nan(...
        length(metrics.time.freq(f).(met).auc.times),...
        length(metrics.time.freq(f).(met).pt),2);  
    metrics.time.freq(f).(met).auc.soz.data = nan(length(metrics.time.freq(f).(met).pt),2);
    % Loop through patients
    for p = 1:length(metrics.time.freq(f).(met).pt)



        % Loop over spike and not a spike
        sp_count = 0;
        for sp = {'spike','not'}
            sp_count = sp_count + 1;
            data = metrics.time.freq(f).(met).pt(p).(sp{1}).data;

            % Get include times
            if strcmp(sp{1},'spike')
                include = include_times(p).times_to_include;
            elseif strcmp(sp{1},'not')
                include = include_times(p).times_to_include_non_spike;
            end

            % Check that data is same size as include times
            if ~isequal(size(data),size(include)), error('what'); end

            % Turn everything not included into a nan
            %if p == 6, error('what'); end
            if rm_rise == 1
                data(~include) = nan;
            end

            % Find first non nan column
            sum_all = sum(data,1);
            non_nan = find(~isnan(sum_all));
            first_non_nan = non_nan(1);
            last_non_nan = non_nan(end);
            first_column = data(:,first_non_nan);
            
            
            metrics.time.freq(f).(met).pt(p).(sp{1}).all_first_last = [data(:,first_non_nan),data(:,last_non_nan)];
            
            % Change all data to be relative to first one
            data = (data-first_column)./abs(first_column);
            last_column = data(:,last_non_nan);

            % AUC for each spike
            if do_cumulative == 1
                auc = nansum(data(:,2:last_non_nan),2); % cumulative sum of relative change
            else
                auc = last_column; % just the last one
            end

            % Take median across all spikes (because a small number are crazy)
            median_auc = nanmedian(auc);
            median_rel_data = nanmedian(data,1);
            iqr_auc = [prctile(auc,25),prctile(auc,75)];
            
            metrics.time.freq(f).(met).pt(p).(sp{1}).all_auc = auc;
            
            % SOZ vs non soz
            if strcmp(sp,'spike') % only makes sense to do it for the spikes
                
                % Get auc for spikes in soz
                auc_soz = auc(logical(is_spike_soz(p).is_soz==1));
                auc_not = auc(logical(is_spike_soz(p).is_soz==0)); % for spikes not in soz
                
                % Take median
                median_auc_soz = median(auc_soz);
                median_auc_not = median(auc_not);
                
                if 0
                    plot(1+randn(length(auc_soz),1)*0.05,auc_soz,'o')
                    hold on
                    plot(2+randn(length(auc_not),1)*0.05,auc_not,'o') 
                end
                
                
                % First versus other
                first = metrics.time.freq(f).(met).pt(p).first;
                other = metrics.time.freq(f).(met).pt(p).other;
                first_auc = (first(:,last_non_nan)-first(:,first_non_nan))./abs(first(:,first_non_nan));
                other_auc = (other(:,last_non_nan)-other(:,first_non_nan))./abs(other(:,first_non_nan));
                first_median = nanmedian(first_auc);
                other_median = nanmedian(other_auc);
                
                metrics.time.freq(f).(met).auc.first_v_other.data(p,:) = [first_median,other_median];
                
                % Depth
                auc_depth = auc(logical(is_spike_depth(p).is_depth));
                auc_not_depth = auc(~logical(is_spike_depth(p).is_depth));
                
                median_auc_depth = median(auc_depth);
                median_auc_not_depth = median(auc_not_depth);
            end
            
            

            % Add it to larger struct
            metrics.time.freq(f).(met).auc.data(p,sp_count) = median_auc;
            metrics.time.freq(f).(met).median_data.data(p,:,sp_count) = median_rel_data;
            metrics.time.freq(f).(met).auc.iqr(p,sp_count,1:2) = iqr_auc;
            
            metrics.time.freq(f).(met).auc.soz.data(p,:) = [median_auc_soz,median_auc_not];
            metrics.time.freq(f).(met).auc.depth.data(p,:) = [median_auc_depth,median_auc_not_depth];
            
            % Pt specific
            metrics.time.freq(f).(met).pt(p).(sp{1}).auc = auc;
            metrics.time.freq(f).(met).pt(p).(sp{1}).temp_data = data;
            
            % Now, to pinpoint the earliest change, loop through times
            for tt = 2:length(metrics.time.freq(f).(met).auc.times)
                if do_cumulative == 1
                    short_auc = nansum(data(:,2:tt),2);
                else
                    short_auc = data(:,tt);
                end
                median_short_auc = nanmedian(short_auc);
                metrics.time.freq(f).(met).auc.short.data(tt,p,sp_count) = median_short_auc;
            end

        end
            
        % Do a two-sample test on the within-patient auc
        sp = metrics.time.freq(f).(met).pt(p).spike.auc;
        not = metrics.time.freq(f).(met).pt(p).not.auc;
        
        if 0
            plot(1+randn(length(sp),1)*0.05,sp,'o')
            hold on
            plot(2+randn(length(sp),1)*0.05,not,'o') 
        end
        
        % t-test
        %[~,p_alt,~,st] = ttest2(sp,not);
        
        % ranksum
        p_alt = ranksum(sp,not);
        
        metrics.time.freq(f).(met).auc.individual_ranksum_p(p) = p_alt; 
        
        % Compare first and last for spikes
        first_and_last = metrics.time.freq(f).(met).pt(p).spike.all_first_last;
        
        if 0
            plot(1+randn(length(first_and_last),1)*0.05,first_and_last(:,1),'o')
            hold on
            plot(2+randn(length(first_and_last),1)*0.05,first_and_last(:,2),'o') 
            pause
            close(gcf)
        end
        
        p_alt = signrank(first_and_last(:,1),first_and_last(:,2));
        metrics.time.freq(f).(met).first_last(p) = p_alt;
    end
    % Do paired ttest on auc
    
    [~,pval,~,st] = ttest(metrics.time.freq(f).(met).auc.data(:,1),...
        metrics.time.freq(f).(met).auc.data(:,2));
    %{
    pval = signrank(metrics.time.freq(f).(met).auc.data(:,1),...
        metrics.time.freq(f).(met).auc.data(:,2));
        %}
    metrics.time.freq(f).(met).auc.pval = pval;
    metrics.time.freq(f).(met).auc.tstat = st.tstat;
    metrics.time.freq(f).(met).auc.short.ps = nan(size(metrics.time.freq(f).(met).auc.short.data,1),1);
    for tt = 1:size(metrics.time.freq(f).(met).auc.short.data,1)
        short_auc = squeeze(metrics.time.freq(f).(met).auc.short.data(tt,:,:));
        [~,p] = ttest(short_auc(:,1),short_auc(:,2));
        metrics.time.freq(f).(met).auc.short.ps(tt) = p;
    end
    
    % Soz
    [~,pval,~,st] = ttest(metrics.time.freq(f).(met).auc.soz.data(:,1),...
        metrics.time.freq(f).(met).auc.soz.data(:,2));
    metrics.time.freq(f).(met).auc.soz.pval = pval;
    metrics.time.freq(f).(met).auc.soz.tstat = st.tstat;
    
    % Depth
    [~,pval,~,st] = ttest(metrics.time.freq(f).(met).auc.depth.data(:,1),...
        metrics.time.freq(f).(met).auc.depth.data(:,2));
    metrics.time.freq(f).(met).auc.depth.pval = pval;
    metrics.time.freq(f).(met).auc.depth.tstat = st.tstat;
    
    % First vs other
    [~,pval,~,st] = ttest(metrics.time.freq(f).(met).auc.first_v_other.data(:,1),...
        metrics.time.freq(f).(met).auc.first_v_other.data(:,2));
    metrics.time.freq(f).(met).auc.first_v_other.pval = pval;
    metrics.time.freq(f).(met).auc.first_v_other.tstat = st.tstat;
end

end