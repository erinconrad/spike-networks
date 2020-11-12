function metrics = remove_early_rise(metrics,pre_spike,wpr,comp_points,rm_rise,alpha)


n_freq_total = 0;
for n = 1:length(metrics)
    n_freq_total = n_freq_total + length(metrics(n).time(1).freq);   
end
adj_alpha = alpha/n_freq_total;

    for n = 1:length(metrics)
        for t = 1:length(metrics(n).time)
            for f = 1:length(metrics(n).time(t).freq)
                curr_freq = metrics(n).time(t).freq(f);
                
                %% Loop through metrics
                fnames = fieldnames(curr_freq);
                
                for fn = 1:length(fnames)
                    if strcmp(fnames{fn},'name'), continue; end
                    met = fnames{fn};
                    curr_met = curr_freq.(met);
                    
                    %% Align times with SD
                    % check that time windows line up
                    if pre_spike(1).windows(t).which ~= metrics(n).time(t).time_window
                        error('Time windows do not align');
                    end
                    
                    if rm_rise == 1
                        
                        % Find where in the SD array the metric array
                        % starts
                        if strcmp(wpr,'cons')
                            ps_windows = pre_spike(1).windows(t).cons_windows;
                        else
                            ps_windows = pre_spike(1).windows(t).all_windows;
                        end
                        shift = find(ps_windows == curr_met.pt(1).times(1));
                        
                        all_t = nan(length(curr_met.pt),length(curr_met.pt(1).times));
                        
                        % Get early rise times
                        if strcmp(wpr,'cons')
                            % Loop over patients
                            for p = 1:length(curr_met.pt)
                                spike_dev = pre_spike(p).windows(t).dev.spike(:,shift:end);
                                
                                % Paired ttest comparing first time window
                                % to subsequent time windows
                                for tt = 2:size(spike_dev,2)
                                    [~,~,~,tstats] = ttest(spike_dev(:,1),spike_dev(:,tt));
                                    all_t(p,tt) = tstats.tstat;
                                    
                                end

                            end
                            
                            % Do one sample test of tstats
                            before_rise_windows = ones(size(all_t,2),1);
                            for tt = 2:size(all_t,2)
                                [~,pval] = ttest(all_t(:,tt));
                                if pval < alpha
                                    before_rise_windows(tt:end) = 0;
                                end
                            end

                        end
                    
                    
                    % Loop over patients
                    for p = 1:length(curr_met.pt)
                        
                        if p > length(curr_met.pt), continue; end

                        curr_pt = curr_met.pt(p);
                        
                        if ~isfield(curr_pt,'spike') || ~isfield(curr_pt,'not')
                            continue;
                        end
                        
                        % Skip it if I don't have pre-spike info for it
                        if p > length(pre_spike), continue; end
                        
                        % check that I have the same pt in pre_spike for
                        % this index
                        if ~strcmp(pre_spike(p).name,curr_pt.name)
                            error('Names do not align');
                        end
                        
                        if strcmp(wpr,'cons')
                            before_rise = repmat(before_rise_windows',...
                                size(curr_pt.spike.data,1),1);
                        else
                            before_rise = pre_spike(p).windows(t).(wpr);
                        end

                        
                        

                      %  if strcmp(met,'sd'), error('look'); end
                        % Reduce spike data to only those times before rise
                        if size(before_rise,1) > size(curr_pt.spike.data,1)
                            before_rise = before_rise(1:size(curr_pt.spike.data,1),:);
                        end
                        
                        if size(before_rise,2) > size(curr_pt.spike.data,2)
                            before_rise = before_rise(:,size(before_rise,2) - size(curr_pt.spike.data,2) +1: end);
                        end
                        curr_pt.spike.data(before_rise==0) = nan;
                        
                        % Get the mode across all spikes (this is what I
                        % will use to reduce the not a spike data)
                        before_rise_mode = mode(before_rise,1);
                        before_rise_mode = repmat(before_rise_mode,...
                            size(curr_pt.not.data,1),1);

                        % Reduce not spike data to only those times before
                        % mode rise
                        curr_pt.not.data(before_rise_mode==0) = nan;
                        curr_met.pt(p) = curr_pt;
                        
                        
                    end
                    end
                        
                        
                        % Loop over spike and not
                    for p = 1:length(curr_met.pt)
                        curr_pt = curr_met.pt(p);
                        snames = fieldnames(curr_pt);
                        for sn = 1:length(snames)
                            if strcmp(snames{sn},'name') || strcmp(snames{sn},'times') || strcmp(snames{sn},'index_windows')
                                continue; 
                            end
                            sp_or_not = curr_pt.(snames{sn});
                            
                            % initialize slopes (as many as there are
                            % spikes)
                            sp_or_not.slopes = zeros(size(sp_or_not.data,1),1);
                            sp_or_not.zs = [];
                            
                            % Get slopes for each spike
                            for s = 1:size(sp_or_not.data,1)
                                data = sp_or_not.data(s,:)';
                                first_non_nan_data = data(~isnan(data));
                                if isempty(first_non_nan_data), continue; end
                                first_non_nan_data = first_non_nan_data(1);
                                if comp_points == 0
                                    zdat = data';
                                elseif comp_points == 1
                                    zdat = ((data-nanmean(data))./nanstd(data))';
                                elseif comp_points == 2
                                    zdat = ((data-first_non_nan_data)./first_non_nan_data)';
                                elseif comp_points == 3
                                    
                                    zdat = ((data-first_non_nan_data)./nanstd(data))';
                                elseif comp_points == 4
                                    zdat = (data./nanstd(data))';
                                elseif comp_points == 5
                                    zdat = (data-first_non_nan_data)';
                                end
                                sp_or_not.zs = [sp_or_not.zs;zdat];
                                
                                % remove nans
                                zdat(isnan(zdat)) = [];
                                z = zdat';
                                
                                
                                % get z scores
                                %z = (data-mean(data))./std(data);
                                %z = data;
                                
                                % do linear regression
                                x = [ones(length(z),1), (1:length(z))'];
                                y = z;
                                b = x\y;
                                slope = b(2);
                                %sp_or_not.slopes(s) = slope;
                                
                                %sp_or_not.slopes(s) = z(end)-z(1);
                                %if isempty(z), continue; end
                                sp_or_not.slopes(s) = z(end);
                                
                            end
                            sp_or_not.mean_z = nanmean(sp_or_not.zs,1);
                            
                            
                            metrics(n).time(t).freq(f).(met).pt(p).(snames{sn}) = sp_or_not;
                            
                        end
                         
                        
                        
                        % Now compare slopes between spikes and not spikes
                        % with a two sample t-test
                        slopes_sp = metrics(n).time(t).freq(f).(met).pt(p).spike.slopes;
                        slopes_not = metrics(n).time(t).freq(f).(met).pt(p).not.slopes;
                        
                        [~,pval,~,stats] = ttest2(slopes_sp,slopes_not);
                        metrics(n).time(t).freq(f).(met).pt(p).test.p = pval;
                        metrics(n).time(t).freq(f).(met).pt(p).test.stats = stats;
                        metrics(n).time(t).freq(f).(met).pt(p).test.alpha = adj_alpha;
                        metrics(n).time(t).freq(f).(met).pt(p).test.h = p<adj_alpha;
                        metrics(n).time(t).freq(f).(met).pt(p).test.slopes{1} = slopes_sp;
                        metrics(n).time(t).freq(f).(met).pt(p).test.slopes{2} = slopes_not;
                       
                       
                        
                    end
                end
            end
        end
    end
end