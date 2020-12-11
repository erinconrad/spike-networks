function metrics = generate_summary_stats(metrics,met,include_times)

nfreq = length(metrics.time.freq);

for f = 1:nfreq
    
    
    % Prep array
    metrics.time.freq(f).(met).auc.data = zeros(...
        length(metrics.time.freq(f).(met).pt),2);  
    metrics.time.freq(f).(met).auc.info = {'spike','not'};
    metrics.time.freq(f).(met).median_data.data = zeros(...
        length(metrics.time.freq(f).(met).pt),...
        size(metrics.time.freq(f).(met).pt(1).spike.data,2),2); 
    metrics.time.freq(f).(met).median_data.info = {'spike','not'};

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
            data(~include) = nan;

            % Find first non nan column
            sum_all = sum(data,1);
            non_nan = find(~isnan(sum_all));
            first_non_nan = non_nan(1);
            first_column = data(:,first_non_nan );

            % Change all data to be relative to first one
            data = (data-first_column)./abs(first_column);

            % AUC for each spike
            auc = nansum(data,2);

            % Take median across all spikes (because a small number are crazy)
            median_auc = median(auc);
            median_rel_data = median(data,1);

            % Add it to larger struct
            metrics.time.freq(f).(met).auc.data(p,sp_count) = median_auc;
            metrics.time.freq(f).(met).median_data.data(p,:,sp_count) = median_rel_data;


        end
            
        
    end
    % Do paired ttest on auc
    [~,pval,~,st] = ttest(metrics.time.freq(f).(met).auc.data(:,1),...
        metrics.time.freq(f).(met).auc.data(:,2));
    metrics.time.freq(f).(met).auc.pval = pval;
    metrics.time.freq(f).(met).auc.tstat = st.tstat;
        
    
    
end

end