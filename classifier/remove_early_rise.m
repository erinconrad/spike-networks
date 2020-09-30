function metrics = remove_early_rise(metrics,pre_spike)
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
                    
                    % Loop over patients
                    for p = 1:length(curr_met.pt)
                        curr_pt = curr_met.pt(p);
                        
                        % Skip it if I don't have pre-spike info for it
                        if p > length(pre_spike), continue; end
                        
                        % check that I have the same pt in pre_spike for
                        % this index
                        if ~strcmp(pre_spike(p).name,curr_pt.name)
                            error('Names do not align');
                        end
                        
                        % check that time windows line up
                        if pre_spike(p).windows(t).which ~= metrics(n).time(t).time_window
                            error('Time windows do not align');
                        end
                        
                        % Get the time windows before the early spike rise
                        before_rise = pre_spike(p).windows(t).before_rise;
                        
                        % Get the mode across all spikes (this is what I
                        % will use to reduce the not a spike data)
                        before_rise_mode = mode(before_rise,1);
                        before_rise_mode = repmat(before_rise_mode,...
                            size(curr_pt.not.data,1),1);
                        
                        % Reduce spike data to only those times before rise
                        curr_pt.spike.data(before_rise==0) = nan;
                        
                        % Reduce not spike data to only those times before
                        % mode rise
                       % curr_pt.not.data(before_rise_mode==0) = nan;
                        
                        % Loop over spike and not
                        snames = fieldnames(curr_pt);
                        for sn = 1:length(snames)
                            if strcmp(snames{sn},'name'), continue; end
                            sp_or_not = curr_pt.(snames{sn});
                            
                            % initialize slopes (as many as there are
                            % spikes)
                            sp_or_not.slopes = zeros(size(sp_or_not.data,1),1);
                            
                            % Get slopes for each spike
                            for s = 1:size(sp_or_not.data,1)
                                data = sp_or_not.data(s,:)';
                                
                                % remove nans
                                data(isnan(data)) = [];
                                
                                % get z scores
                                z = (data-mean(data))./std(data);
                                %z = data;
                                
                                % do linear regression
                                x = [ones(length(z),1), (1:length(z))'];
                                y = z;
                                b = x\y;
                                slope = b(2);
                                sp_or_not.slopes(s) = slope;
                                
                                %sp_or_not.slopes(s) = z(end)-z(1);
                            end
                            
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