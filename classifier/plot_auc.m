function plot_auc(metrics,met)

jitter_amount = 0.05;

if strcmp(met,'sd')
    nfreq = 1;
else
    nfreq = length(metrics.time.freq);
end

%% Initialize figure
figure
set(gcf,'position',[1 100 1450 270])
[ha, pos] = tight_subplot(1, nfreq, [0.10 0.05], [0.12 0.11], [0.11 0.01]);


for f = 1:nfreq
    axes(ha(f))
    
    %% Get data
    dat_sp = metrics.time.freq(f).(met).auc.data(:,1);
    dat_not = metrics.time.freq(f).(met).auc.data(:,2);
    iqr_sp = squeeze(metrics.time.freq(f).(met).auc.iqr(:,1,:));
    iqr_not = squeeze(metrics.time.freq(f).(met).auc.iqr(:,2,:));
    pval = metrics.time.freq(f).(met).auc.pval;
    np = length(dat_sp);
    
    %% Do the plot
    for p = 1:np
        
        % Plot spike
        errorbar(p-0.1,dat_sp(p),...
            dat_sp(p) - iqr_sp(p,1),iqr_sp(p,2) - dat_sp(p),'ro',...
            'markersize',10);
        hold on
        
        % Plot not spike
        % Plot spike
        errorbar(p+0.1,dat_not(p),...
            iqr_not(p,2) - dat_not(p),dat_not(p) - iqr_not(p,1),'ko',...
            'markersize',10);
        
        

    end
    

    %% Formatting
    if f == 2
        xlabel('Patient');
    end
    xticklabels([])
    if f == 2
        legend('Spike','No spike','fontsize',20,'location','southwest')
    end
    if f == 1
        ylabel(sprintf('Relative %s summed\nacross pre-spike epoch',met));
    end
    title(sprintf(metrics.time.freq(f).name))
    set(gca,'fontsize',20)
    
    %% Add p-value
    yl = get(gca,'ylim');
    xl = get(gca,'xlim');
    p_xloc = 0.5*(xl(1)+xl(2));
    p_yloc = yl(1) + 0.9*(yl(2)-yl(1));
    pp = pretty_p(pval,3);
    text(p_xloc,p_yloc,pp,'horizontalalignment','center','fontsize',20);
    

    
end



end

function out = add_jitter(sz,jitter_amount)
    out = randn(sz,1)*jitter_amount;

end