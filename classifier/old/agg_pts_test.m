function metrics = agg_pts_test(metrics)

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
                
                % Initialize arrays for all
                % patients
                all_z_spike = [];
                all_z_not = [];
                all_auc_spike = [];
                all_auc_not= [];
                
                % Loop over patients
                for p = 1:length(curr_met.pt)
                    
                    if ~isfield(curr_met.pt(p),'test'), continue; end

                    
                    %all_z_spike = [all_z_spike;curr_met.pt(p).spike.mean_z];
                    %all_z_not = [all_z_not;curr_met.pt(p).not.mean_z];
                    
                    all_z_spike = [all_z_spike;curr_met.pt(p).spike.median_z];
                    all_z_not = [all_z_not;curr_met.pt(p).not.median_z];
                    
                    all_auc_spike = [all_auc_spike;median(curr_met.pt(p).spike.auc)];
                    all_auc_not = [all_auc_not;median(curr_met.pt(p).not.auc)];
                   % if strcmp(met,'sd'), error('look'); end
                end
                
                
                
                curr_met.all_z_spike = all_z_spike;
                curr_met.all_z_not = all_z_not;
                
                curr_met.auc_spike = all_auc_spike;
                curr_met.auc_not = all_auc_not;
                
                [~,curr_met.auc_pval] = ttest(all_auc_spike,all_auc_not);
                metrics(n).time(t).freq(f).(met) = curr_met;
            end
            
        end
        
    end
    
end


end