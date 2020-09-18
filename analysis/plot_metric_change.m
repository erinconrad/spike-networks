function plot_metric_change(stats,windows,out_folder,paired,is_spike,metric)

network_count = length(stats);
time_count = length(stats(1).time);
n_freq_abs = length(stats(1).time(1).freq);

figure
set(gcf,'position',[1 100 1500 length(windows)*200+100])
[ha, ~] = tight_subplot(length(windows), n_freq_abs+1, [0.08 0.02], [0.2 0.12], [0.06 0.005]);
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
            
            % Get z scores and h's
            z_curr = stats(n).time(t).freq(f).(metric).z_curr;
            h = stats(n).time(t).freq(f).(metric).tests.(paired).h(is_spike);
            np = size(z_curr,1);
            
            % Loop over patients and plot
            for i = 1:np
                plot(times,squeeze(z_curr(i,:,is_spike)),'ko'); 
                hold on
            end
            
            % Get overall trend line
            x = repmat(times',size(z_curr,1),1);
            y = z_curr(:,:,is_spike);
            y = y(:);
            x = x(:);
            X = [ones(length(x),1),x];
            b = X\y;
            
            if h == 1
                plot(times',b(1)+b(2)*times,'g','linewidth',3);
            else
                plot(times',b(1)+b(2)*times,'k','linewidth',3);
            end
            
            set(gca,'fontsize',20)
            if f == 4 && tcount == length(windows)
                xlabel('Time relative to spike peak (s)')
            end 
            if n == 2
                if length(windows) == 3 && tcount == 2
                    ylabel(sprintf('metric (z-score)'))
                elseif length(windows) == 1
                    ylabel(sprintf('metric (z-score)'))
                end
            end
            
            % adjust z_range
            if max(max(z_curr(:,:,1))) > z_range(tcount,2)
                z_range(tcount,2) = max(max(z_curr(:,:,1)));
            end

            if min(min(z_curr(:,:,1))) < z_range(tcount,1)
                z_range(tcount,1) = min(min(z_curr(:,:,1)));
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
    t = ceil(sp/(n_freq_abs+1));
    ylim(z_range(t,:))
    
    if mod(sp,(n_freq_abs+1)) ~= 1
        yticklabels([])
    end
end

print(gcf,[out_folder,'updated_net_change_',num2str(is_spike)],'-depsc');