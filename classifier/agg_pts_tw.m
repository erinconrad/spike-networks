function agg_pts_tw(stats,met,windows,method)

network_count = length(stats);
time_count = length(stats(1).time);
n_freq_abs = length(stats(1).time(1).freq);

figure
set(gcf,'position',[1 100 1500 length(windows)*200+100])
[ha, ~] = tight_subplot(length(windows), n_freq_abs+1, [0.08 0.01], [0.2 0.12], [0.07 0.005]);
z_range = zeros(length(windows),2);

for n = 1:network_count

    net_name = stats(n).name;
    tcount = 0;
    
    for t = 1:time_count
        
        % Skip it if I want
        if ~ismember(stats(n).time(t).time_window,windows), continue; end
        tcount = tcount+1;
        times = stats(n).time(t).times;
        
        nfreq = length(stats(n).time(t).freq);
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
            sp = (n_freq_abs+1)*(tcount-1) + f + column_add;
            axes(ha(sp));
            
            % Get stats
            dat_sp = stats(n).time(t).freq(f).(met).all_z_spike;
            dat_not = stats(n).time(t).freq(f).(met).all_z_not;
            %h = stats(n).time(t).freq(f).(met).pt(pt_id).(method).h;
            pval = stats(n).time(t).freq(f).(met).(method).p;
            text_p = get_asterisks(pval,(n_freq_abs+1));
            
            dat_sp_mean = nanmean(dat_sp,1);
            dat_not_mean = nanmean(dat_not,1);
            dat_sp_se = nanstd(dat_sp,0,1)./sqrt(sum(~isnan(dat_sp),1));
            dat_not_se = nanstd(dat_not,0,1)./sqrt(sum(~isnan(dat_not),1));
            
            errorbar(times,dat_sp_mean,dat_sp_se);
            hold on
            errorbar(times,dat_not_mean,dat_not_se);
            all_slopes = [dat_sp(:);dat_not(:)];
            yl = get(gca,'ylim');
            line_height = max(all_slopes) + yl(2)+0.2*(yl(2)-yl(1));
           
            set(gca,'fontsize',20)
            
            
            % adjust z_range
            if max(all_slopes) > z_range(tcount,2)
                z_range(tcount,2) = max(all_slopes);
            end

            if min(all_slopes) < z_range(tcount,1)
                z_range(tcount,1) = min(all_slopes);
            end
            
            if n == 2
                if length(windows) == 3 && tcount == 2
                    ylabel(sprintf('%s slope\n(z-score)',met))
                elseif length(windows) == 1
                    ylabel(sprintf('%s slope\n(z-score)',met))
                end
            end
            
            if tcount == 1 && strcmp(net_name,'coherence') == 1
                title(sprintf('%s',...
                    strrep(stats(n).time(t).freq(f).name,'_',' ')))
            elseif tcount == 1 && strcmp(net_name,'simple') == 1
                title('correlation')
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

