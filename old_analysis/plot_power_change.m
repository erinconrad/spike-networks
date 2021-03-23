function plot_power_change(sig_dev,windows,out_folder,paired,is_spike)

time_count = length(sig_dev);

figure
set(gcf,'position',[1 100 1399 length(windows)*200+100])
[ha, ~] = tight_subplot(length(windows), 1, [0.1 0.01], [0.2 0.1], [0.06 0.02]);

tcount = 0;
for t = 1:time_count
    
    % Skip it if I want
    if ~ismember(sig_dev(t).is_spike(1).time_window,windows), continue; end
    tcount = tcount+1;
    axes(ha(tcount));
    
    s = is_spike; 
    
    z_curr = sig_dev(t).is_spike(s).z_curr;
    times = sig_dev(t).is_spike(s).times;
    h = sig_dev(t).tests.(paired).h(is_spike);
    np = size(z_curr,1);
    % Loop over patients and plot
    for i = 1:np
        plot(times,squeeze(z_curr(i,:)),'ko');
        hold on
    end
    
    % Plot an overall trend line
    x = repmat(times',size(z_curr,1),1);
    y = z_curr(:);
    x = x(:);
    x = [ones(length(x),1),x];
    b = x\y;
    
    if h == 1
        plot(times,b(1)+b(2)*times,'g','linewidth',3);
    else
        plot(times,b(1)+b(2)*times,'k','linewidth',3);
    end
    
    set(gca,'fontsize',20)
    
    xlabel('Time relative to spike peak (s)')
    if length(windows) == 3 && tcount == 2
        ylabel(sprintf('Power change from\nfirst time (z-score)'))
    elseif length(windows) == 1
        ylabel(sprintf('Power change from\nfirst time (z-score)'))
    end
    
end


print(gcf,[out_folder,'updated_power_change_',num2str(is_spike)],'-depsc');
end