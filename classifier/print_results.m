function print_results(metrics,met)

if contains(met,'sd')
    nfreq = 1;
else
    nfreq = 3;
end

if contains(met,'sd')
    pretty_name = 'power';
elseif contains(met,'ers')
    pretty_name = 'power';
elseif strcmp(met,'ns_avg')
    pretty_name = 'node strength averaged across electrodes';
elseif strcmp(met,'ns_auto')
    pretty_name = 'node strength in the peak IED electrode';
end

for f = 1:nfreq
    
    if contains(met,'sd')
        freq_name = '';
    else
        freq_name = [metrics.time.freq(f).name];
    end
    
    if metrics.time.freq(f).(met).auc.pval < 0.05/nfreq
        fprintf(['\nThe pre-IED relative %s %s change was higher for IED periods (M = %1.3f, SD = %1.3f) than IED-free periods (M = %1.3f, SD = %1.3f) (t(%d) = %1.1f, p = %1.3f)\n'],...
            freq_name,pretty_name,...
            mean(metrics.time.freq(f).(met).auc.data(:,1)),std(metrics.time.freq(f).(met).auc.data(:,1)),...
            mean(metrics.time.freq(f).(met).auc.data(:,2)),std(metrics.time.freq(f).(met).auc.data(:,2)),...
            metrics.time.freq(f).(met).auc.df,metrics.time.freq(f).(met).auc.tstat,...
            metrics.time.freq(f).(met).auc.pval);
    else
        fprintf(['\nThere was no difference in the pre-IED relative %s %s change in IED periods (M = %1.3f, SD = %1.3f) compared to IED-free periods (M = %1.3f, SD = %1.3f) (t(%d) = %1.1f, p = %1.3f)\n'],...
            freq_name,pretty_name,...
            mean(metrics.time.freq(f).(met).auc.data(:,1)),std(metrics.time.freq(f).(met).auc.data(:,1)),...
            mean(metrics.time.freq(f).(met).auc.data(:,2)),std(metrics.time.freq(f).(met).auc.data(:,2)),...
            metrics.time.freq(f).(met).auc.df,metrics.time.freq(f).(met).auc.tstat,...
            metrics.time.freq(f).(met).auc.pval);
        
    end
    
    if metrics.time.freq(f).(met).auc.pval < 0.05/nfreq
        fprintf('\nThe earliest significant pre-IED relative %s %s change occurred %d ms prior to the IED peak\n',...
            freq_name,pretty_name,...
            -metrics.time.freq(f).(met).auc.short.first_sig_time*1000);
    end
    
    
end


end