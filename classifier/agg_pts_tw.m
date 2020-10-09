function agg_pts_tw(stats,met,windows,method)

if strcmp(met,'sd')
    met_text = 'power';
else
    met_text = met;
end
network_count = length(stats);
time_count = length(stats(1).time);
n_freq_abs = length(stats(1).time(1).freq);
if strcmp(met,'sd'), n_freq_abs = 0; end

figure
if strcmp(met,'sd')
    set(gcf,'position',[1 100 600 length(windows)*200+200])
    [ha, ~] = tight_subplot(length(windows), 1, [0.08 0.01], [0.2 0.12], [0.10 0.09]);
else
    set(gcf,'position',[1 100 1500 length(windows)*200+100])
    [ha, ~] = tight_subplot(length(windows), n_freq_abs+1, [0.08 0.01], [0.2 0.12], [0.05 0.005]);
end
z_range = zeros(length(windows),2);

if strcmp(met,'sd')
    network_count = 1;
end
for n = 1:network_count

    net_name = stats(n).name;
    tcount = 0;
    
    for t = 1:time_count
        
        % Skip it if I want
        if ~ismember(stats(n).time(t).time_window,windows), continue; end
        tcount = tcount+1;
        times = stats(n).time(t).freq(1).(met).pt(1).times;
        
        nfreq = length(stats(n).time(t).freq);
        if strcmp(met,'sd')
            nfreq = 1;
        end
        for f = 1:nfreq
            
            % Get appropriate subplot
            if strcmp(net_name,'coherence') == 1
                column_add = 1;
            else
                column_add = 0;
            end
            % this adds the number of frequencies + 1 if it's on the 2nd
            % time point (to move down a row), and it adds which frequency
            % (which is 1 if simple) and adds 1 if coherence, to start with
            % the 2nd column for coherence
            %sp = f + column_add;
            if strcmp(met,'sd')
                sp = tcount;
            else
                sp = (n_freq_abs+1)*(tcount-1) + f + column_add;
            end
            axes(ha(sp));
            
            % Get stats
            dat_sp = stats(n).time(t).freq(f).(met).all_z_spike;
            dat_not = stats(n).time(t).freq(f).(met).all_z_not;
            %h = stats(n).time(t).freq(f).(met).pt(pt_id).(method).h;
            pval = stats(n).time(t).freq(f).(met).(method).p;
            text_p = get_asterisks(pval,(n_freq_abs+1));
            
            % Remove columns with only one non nan
            for j = 1:size(dat_sp,2)
                if sum(~isnan(dat_sp(:,j))) == 1
                    dat_sp(:,j) = nan;
                end
                
                if sum(~isnan(dat_not(:,j))) == 1
                    dat_not(:,j) = nan;
                end
            end
            
            dat_sp_mean = nanmean(dat_sp,1);
            dat_not_mean = nanmean(dat_not,1);
            dat_sp_se = nanstd(dat_sp,0,1)./sqrt(sum(~isnan(dat_sp),1));
            dat_not_se = nanstd(dat_not,0,1)./sqrt(sum(~isnan(dat_not),1));
            
            
            
            errorbar(times,dat_sp_mean,dat_sp_se,'linewidth',2);
            hold on
            errorbar(times,dat_not_mean,dat_not_se,'linewidth',2);
            all_slopes = [dat_sp(:);dat_not(:)];
            
            %% Significance testing 
            % Test just the last point
            all_nan_columns = sum(isnan(dat_sp),1) == size(dat_sp,1);
            last_non_nan = find(all_nan_columns);
            last_non_nan(last_non_nan == 1) = [];
            last_non_nan = last_non_nan(1)-1;
            [~,pval] = ttest(dat_sp(:,last_non_nan),dat_not(:,last_non_nan));
            prettyp = pretty_p(pval,n_freq_abs+1);
            yl = get(gca,'ylim');
            xl = get(gca,'xlim');
            
            maxyloc = max([dat_sp_mean(last_non_nan)+dat_sp_se(last_non_nan),...
                dat_not_mean(last_non_nan)+dat_not_se(last_non_nan)]);
                
            line_height = maxyloc + 0.1*(yl(2)-yl(1));
            text((xl(1)+xl(2))/2,line_height,...
                sprintf('%s',prettyp),'fontsize',20,...
                'horizontalalignment','center')
            ylim([yl(1) line_height + 0.1*(yl(2)-yl(1))]);
            
            
            
           
            set(gca,'fontsize',20)
            
            
            % adjust z_range
            if max(all_slopes) > z_range(tcount,2)
                z_range(tcount,2) = max(all_slopes);
            end

            if min(all_slopes) < z_range(tcount,1)
                z_range(tcount,1) = min(all_slopes);
            end
            
            if n == network_count
                if length(windows) == 3 && tcount == 2
                    ylabel(sprintf('Normalized %s',met_text))
                elseif length(windows) == 1
                    ylabel(sprintf('Normalized %s',met_text))
                end
            end
            
            if ~strcmp(met,'sd')
                if tcount == 1 && strcmp(net_name,'coherence') == 1
                    title(sprintf('%s',...
                        strrep(stats(n).time(t).freq(f).name,'_',' ')))
                elseif tcount == 1 && strcmp(net_name,'simple') == 1
                    title('correlation')
                end
            end
            
            if tcount == length(windows)
                if strcmp(met,'sd') || (f == floor((n_freq_abs+1)/2) && n == 1) 
                    xlabel('Time relative to spike peak (s)')
                end
            end
        end
        
        
    end
    
    
end

for sp = 1:length(ha)
    axes(ha(sp))
    %t = ceil(sp/(n_freq_abs+1));
    %ylim(z_range(t,:))
    
    if mod(sp,(n_freq_abs+1)) ~= 1
        yticklabels([])
    end
end

end

