function metrics = soz_comparison(metrics,met,out_folder)

jitter_amount = 0.05;

%% Pretty names
if strcmp(met,'sd')
    pretty_name = 'absolute\newlinepower';
else
    pretty_name = met;
end

figure

if strcmp(met,'sd')
    nfreq = 1;
    set(gcf,'position',[1 100 400 270])
    [ha, pos] = tight_subplot(1, nfreq, [0.10 0.05], [0.12 0.11], [0.22 0.01]);
else
    nfreq = length(metrics.time.freq);
    set(gcf,'position',[1 100 1100 270])
    [ha, pos] = tight_subplot(1, nfreq, [0.10 0.05], [0.12 0.11], [0.11 0.01]);
end



for f = 1:nfreq
    axes(ha(f))
    dat_soz = metrics.time.freq(f).(met).auc.soz.data(:,1);
    dat_not = metrics.time.freq(f).(met).auc.soz.data(:,2);
    pval = metrics.time.freq(f).(met).auc.soz.pval;
    np = length(dat_soz);
    x_soz = ones(np,1) + add_jitter(np,jitter_amount);
    x_not = 2*ones(np,1) + add_jitter(np,jitter_amount);
    
    plot(x_soz,dat_soz,'ro','markersize',10,'linewidth',2)
    hold on
    plot(x_not,dat_not,'ko','markersize',10,'linewidth',2)
    
    xticks([1 2])
    xticklabels({'SOZ','Not SOZ'})
    if f == 1
        ylabel(sprintf('Pre-IED %s change',pretty_name))
    end
    xlim([0.5 2.5])
    set(gca,'fontsize',20)
    
    p_pretty = pretty_p(pval,nfreq);
    yl = get(gca,'ylim');
    p_xloc = 1.5;
    p_yloc = yl(1) + 0.9*(yl(2)-yl(1));
    text(p_xloc,p_yloc,p_pretty,'HorizontalAlignment','Center','fontsize',20)
    
    
    
end

print(gcf,[out_folder,sprintf('soz_%s',met)],'-dpng')

end

function out = add_jitter(n,amount)
    out = randn(n,1)*amount;

end