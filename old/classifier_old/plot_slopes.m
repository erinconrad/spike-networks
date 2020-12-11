function plot_slopes(stats,met,pt_id,windows)

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
            
            % Get slopes
            slopes = stats(n).time(t).freq(f).(met).pt(pt_id).test.slopes;
            h = stats(n).time(t).freq(f).(met).pt(pt_id).test.h;
            pval = stats(n).time(t).freq(f).(met).pt(pt_id).test.p;
            alpha = stats(n).time(t).freq(f).(met).pt(pt_id).test.alpha;
            text_p = get_asterisks(pval,(n_freq_abs+1));
            
            plot(1,slopes{1},'ko');
            hold on
            plot(2,slopes{2},'ko');
            all_slopes = [slopes{1};slopes{2}];
            yl = get(gca,'ylim');
            line_height = max(all_slopes) + yl(2)+0.2*(yl(2)-yl(1));
            plot([1 2],[line_height line_height],...
                'k');
            text(1.5,line_height+0.2*(yl(2)-yl(1)),sprintf('p = %1.3f%s',pval,text_p),'HorizontalAlignment','Center',...
                'Fontsize',20);
            xlim([0.3 2.7])
            yl = get(gca,'ylim');
            yl(2) = line_height+0.2*(yl(2)-yl(1));
            ylim(yl);
            xticks([1 2])
            xticklabels({'Spike','Not'})
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