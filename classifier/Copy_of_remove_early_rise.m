function metrics = remove_early_rise(metrics,pre_spike,wpr,comp_points,rm_rise)
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
                
                %{
                %% Reduce sp diff
                sp_diff = curr_freq.sp_diff;
                
                % Loop over pts
                for p = 1:length(sp_diff.F.pt)
                    
                    % Get the time windows before the early spike rise
                    before_rise = pre_spike(p).windows(t).before_rise;
                    
                    for fn = {'F','score'}
                        sp_sub = sp_diff.(fn{1}).pt(p);
                        sp_sub.spike.data(before_rise == 0) = nan;
                        
                        % Get slopes
                        for sp = {'spike','not'}
                            sp_sub_2 = sp_sub.(sp{1});
                            data = sp_sub_2.data;
                            % DO THIS*******
                            sp_sub_2.slopes = zeros(size(data,1),1);
                            
                            % Loop over spikes
                            for s = 1:size(data,1)
                                
                                sdata = squeeze(data(s,:));
                                
                                % remove nans
                                sdata(isnan(sdata)) = [];
                                
                                z = (sdata-mean(sdata))./std(sdata);
                                %z = data;
                                
                                % do linear regression
                                x = [ones(length(z),1), (1:length(z))'];
                                y = z';
                                b = x\y;
                                slope = b(2);
                                sp_sub_2.slopes(s) = slope;
                                
                                
                            end
                            
                            metrics(n).time(t).freq(f).sp_diff.(fn{1}).pt(p).(sp{1}).slopes = sp_sub_2.slopes;
                        end
                        
                        % Now compare slopes in sp vs not using two sample
                        % ttest
                        slopes_sp = metrics(n).time(t).freq(f).sp_diff.(fn{1}).pt(p).spike.slopes;
                        slopes_not = metrics(n).time(t).freq(f).sp_diff.(fn{1}).pt(p).not.slopes;
                        [~,pval,~,stats1] = ttest2(slopes_sp,slopes_not);
                        metrics(n).time(t).freq(f).sp_diff.(fn{1}).pt(p).test.p = pval;
                        metrics(n).time(t).freq(f).sp_diff.(fn{1}).pt(p).test.stats = stats1;
                        metrics(n).time(t).freq(f).sp_diff.(fn{1}).pt(p).test.slopes{1} = slopes_sp;
                        metrics(n).time(t).freq(f).sp_diff.(fn{1}).pt(p).test.slopes{2} = slopes_not;
                    end

                    if 0
                        figure
                        subplot(2,2,1)
                        plot(curr_freq.sp_diff.F.pt(p).spike.data','ko')
                        title('Spike F')
                        
                        subplot(2,2,2)
                        plot(curr_freq.sp_diff.F.pt(p).not.data','ko')
                        title('Not spike F')
                        
                        subplot(2,2,3)
                        plot(curr_freq.sp_diff.score.pt(p).spike.data','ko')
                        title('Spike score')
                        
                        subplot(2,2,4)
                        plot(curr_freq.sp_diff.score.pt(p).not.data','ko')
                        title('Not spike score')
                        pause
                        close gcf
                    end
                end
                %}
                
                %% Loop through metrics
                fnames = fieldnames(curr_freq);
                
                for fn = 1:length(fnames)
                    if strcmp(fnames{fn},'name'), continue; end
                    met = fnames{fn};
                    curr_met = curr_freq.(met);
                    
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
                        
                        % check that time windows line up
                        if pre_spike(p).windows(t).which ~= metrics(n).time(t).time_window
                            error('Time windows do not align');
                        end
                        
                        %metrics(n).time(t).times = pre_spike(p).windows(t).all_windows;
                        metrics(n).time(t).times = curr_pt.times;
                        

                    %    if p == 10 && strcmp(met,'ers'), error('look'); end
                        
                    
                        if rm_rise == 1
                            
                            % Find where in the SD array the metric array
                            % starts
                            shift = find(pre_spike(p).windows(t).all_windows == curr_pt.times(1));
                            
                            if strcmp(wpr,'cons')
                                % Need to do significance testing on SD
                                % relative to the shifted data
                                spike_dev = pre_spike(p).windows(t).dev.spike(:,shift:end);
                                
                            end
                            
                            if ~isequal(size(before_rise),size(curr_pt.spike.data))
                                fprintf('\nWarning, SD size does not overlap with metric\n');
                                %continue

                                if size(before_rise,2) < size(curr_pt.spike.data,2)
                                    error('what');
                                end
                            end

                            

                             % Get the time windows before the early spike rise
                            before_rise = pre_spike(p).windows(t).(wpr);

                            % Get the mode across all spikes (this is what I
                            % will use to reduce the not a spike data)
                            before_rise_mode = mode(before_rise,1);
                            before_rise_mode = repmat(before_rise_mode,...
                                size(curr_pt.not.data,1),1);

                            % Reduce spike data to only those times before rise
                            curr_pt.spike.data(before_rise(:,shift:end)==0) = nan;

                            % Reduce not spike data to only those times before
                            % mode rise
                            curr_pt.not.data(before_rise_mode(:,shift:end)==0) = nan;
                        end
                        
                        
                        % Loop over spike and not
                        snames = fieldnames(curr_pt);
                        for sn = 1:length(snames)
                            if strcmp(snames{sn},'name'), continue; end
                            sp_or_not = curr_pt.(snames{sn});
                            
                            % initialize slopes (as many as there are
                            % spikes)
                            sp_or_not.slopes = zeros(size(sp_or_not.data,1),1);
                            sp_or_not.zs = [];
                            
                            % Get slopes for each spike
                            for s = 1:size(sp_or_not.data,1)
                                data = sp_or_not.data(s,:)';
                                first_non_nan_data = data(~isnan(data));
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