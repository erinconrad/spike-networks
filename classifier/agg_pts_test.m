function metrics = agg_pts_test(metrics)

alpha = 0.05;

n_freq_total = 0;
for n = 1:length(metrics)
    n_freq_total = n_freq_total + length(metrics(n).time(1).freq);   
end
adj_alpha = alpha/n_freq_total;

for n = 1:length(metrics)
    for t = 1:length(metrics(n).time)
        for f = 1:length(metrics(n).time(t).freq)
            curr_freq = metrics(n).time(t).freq(f);

            % Loop through metrics
            fnames = fieldnames(curr_freq);

            for fn = 1:length(fnames)
                if strcmp(fnames{fn},'name'), continue; end
                met = fnames{fn};
                curr_met = curr_freq.(met);
                
                % Initialize arrays of t statistics and p values for all
                % patients
                all_t = [];
                all_p = [];
                all_slopes = [];
                all_abs_slopes = [];
                all_z_spike = [];
                all_z_not = [];
                
                % Loop over patients
                for p = 1:length(curr_met.pt)
                    
                    if ~isfield(curr_met.pt(p),'test'), continue; end
                    
                    % add t stat and p value to arrays
                    all_t = [all_t;curr_met.pt(p).test.stats.tstat];
                    all_p = [all_p;curr_met.pt(p).test.p];
                    
                    % Add average slopes for both spike and not a spike
                    all_slopes = [all_slopes;...
                        mean((curr_met.pt(p).test.slopes{1})),...
                        mean((curr_met.pt(p).test.slopes{2}))];
                    
                    all_abs_slopes = [all_abs_slopes;...
                        mean(abs(curr_met.pt(p).test.slopes{1})),...
                        mean(abs(curr_met.pt(p).test.slopes{2}))];
                    
                    all_z_spike = [all_z_spike;curr_met.pt(p).spike.mean_z];
                    all_z_not = [all_z_not;curr_met.pt(p).not.mean_z];
                    
                end
                
                % Do 3 methods (one sample ttest, Fisher's method, paired t of mean slopes) to
                % aggregate stats across patients
                fisher_p = fisher_p_value(all_p);
                fisher_h = fisher_p < adj_alpha;
                
                [~,tstat_p,~,tstat_stats] = ttest(all_t);
                tstat_h = tstat_p < adj_alpha;
                
                if isempty(all_slopes), continue; end
                [~,tstatp_p,~,tstatp_stats] = ttest(all_slopes(:,1),all_slopes(:,2));
                tstatp_h = tstatp_p < adj_alpha;
                
                [~,tstatpabs_p,~,tstatpabs_stats] = ttest(all_abs_slopes(:,1),all_abs_slopes(:,2));
                tstatpabs_h = tstatpabs_p < adj_alpha;
                
                % Add to structure
                curr_met.fisher.p = fisher_p;
                curr_met.fisher.h = fisher_h;
                curr_met.fisher.all_p = all_p;
                curr_met.fisher.all_slopes = all_slopes;
                
                curr_met.ttest.p = tstat_p;
                curr_met.ttest.stats = tstat_stats;
                curr_met.ttest.h = tstat_h;
                curr_met.ttest.all_t = all_t;
                curr_met.ttest.all_slopes = all_slopes;
                
                curr_met.ttestp.p = tstatp_p;
                curr_met.ttestp.stats = tstatp_stats;
                curr_met.ttestp.h = tstatp_h;              
                curr_met.ttestp.all_slopes = all_slopes;
                
                curr_met.ttestpabs.p = tstatpabs_p;
                curr_met.ttestpabs.stats = tstatpabs_stats;
                curr_met.ttestpabs.h = tstatpabs_h;              
                curr_met.ttestpabs.all_slopes = all_abs_slopes;
                
                curr_met.all_z_spike = all_z_spike;
                curr_met.all_z_not = all_z_not;
                
                metrics(n).time(t).freq(f).(met) = curr_met;
            end
            
        end
        
    end
    
end


end