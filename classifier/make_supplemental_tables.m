function make_supplemental_tables(meta_metrics,results_folder,do_save)

metrics_sd = meta_metrics.metrics_sd;
metrics_ers = meta_metrics.metrics_ers;
metrics_ns_big= meta_metrics.metrics_ns_big;
metrics_ns_avg = meta_metrics.metrics_ns_avg;

which_mets = {'sd_auto','ers_auto','ns_auto','ns_avg'};
all_metrics = [metrics_sd,metrics_ers,metrics_ns_big,metrics_ns_avg];
pretty_met_names = {'absolute power','power','peak-IED electrode node strength','electrode average node strength'};
nfreq = [1 3 3 3];

%% Make a table for SOZ
soz_means = {};
non_soz_means = {};
soz_t = {};
soz_p = {};
all_names = {};

% Print SOZ stuff
fprintf('\n');
for m = 1:length(which_mets)
    met = which_mets{m};
    %fprintf('\n%s\n',met);
    
    metrics = all_metrics(m);
    nf = nfreq(m);
    for f = 1:nf
        if contains(met,'sd')
            fname = '';
        else
            fname = [metrics.time.freq(f).name,' '];
        end
        
        all_names = [all_names;[fname,pretty_met_names{m}]];
        soz_means = [soz_means;sprintf('%1.2f (%1.2f)',...
            nanmean(metrics.time.freq(f).(met).auc.soz.data(:,1)),...
            nanstd(metrics.time.freq(f).(met).auc.soz.data(:,1)))];
        
        non_soz_means = [non_soz_means;sprintf('%1.2f (%1.2f)',...
            nanmean(metrics.time.freq(f).(met).auc.soz.data(:,2)),...
            nanstd(metrics.time.freq(f).(met).auc.soz.data(:,2)))];
        
        soz_t = [soz_t;sprintf('t(%d) = %1.1f',...
            metrics.time.freq(f).(met).auc.soz.df,...
            metrics.time.freq(f).(met).auc.soz.tstat)];
        
        soz_p = [soz_p;sprintf('p = %1.3f',...
            metrics.time.freq(f).(met).auc.soz.pval)];
        
    end
end

soz_table = table(soz_means,non_soz_means,soz_t,soz_p);
if do_save
writetable(soz_table,[results_folder,'tables/','soz_table.csv']);
end

%% Make a table for lead vs non-lead
lead_means = {};
non_lead_means = {};
lead_t = {};
lead_p = {};
all_names = {};

% Print lead stuff
fprintf('\n');
for m = 1:length(which_mets)
    met = which_mets{m};
    %fprintf('\n%s\n',met);
    
    metrics = all_metrics(m);
    nf = nfreq(m);
    for f = 1:nf
        if contains(met,'sd')
            fname = '';
        else
            fname = [metrics.time.freq(f).name,' '];
        end
        
        all_names = [all_names;[fname,pretty_met_names{m}]];
        lead_means = [lead_means;sprintf('%1.2f (%1.2f)',...
            nanmean(metrics.time.freq(f).(met).auc.first_v_other.data(:,1)),...
            nanstd(metrics.time.freq(f).(met).auc.first_v_other.data(:,1)))];
        
        non_lead_means = [non_lead_means;sprintf('%1.2f (%1.2f)',...
            nanmean(metrics.time.freq(f).(met).auc.first_v_other.data(:,2)),...
            nanstd(metrics.time.freq(f).(met).auc.first_v_other.data(:,2)))];
        
        lead_t = [lead_t;sprintf('t(%d) = %1.1f',...
            metrics.time.freq(f).(met).auc.first_v_other.df,...
            metrics.time.freq(f).(met).auc.first_v_other.tstat)];
        
        lead_p = [lead_p;sprintf('p = %1.3f',...
            metrics.time.freq(f).(met).auc.first_v_other.pval)];
        
    end
end


lead_table = table(lead_means,non_lead_means,lead_t,lead_p);
if do_save
    writetable(lead_table,[results_folder,'tables/','lead_table.csv']);
end



end