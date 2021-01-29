function plot_auc(metrics,met,out_folder)

do_all_pts = 0;
jitter_amount = 0.05;

%% Pretty names
if strcmp(met,'sd')
    pretty_name = 'absolute\newlinepower';
else
    pretty_name = met;
end

if do_all_pts

    

    %% Initialize figure
    figure
    if strcmp(met,'sd')
        nfreq = 1;
        set(gcf,'position',[1 100 650 270])
        [ha, pos] = tight_subplot(1, nfreq, [0.10 0.05], [0.12 0.11], [0.15 0.01]);
    else
        nfreq = length(metrics.time.freq);
        set(gcf,'position',[1 100 1450 270])
        [ha, pos] = tight_subplot(1, nfreq, [0.10 0.05], [0.12 0.11], [0.11 0.01]);
    end




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
                'markersize',15,'linewidth',2);
            hold on

            % Plot not spike
            % Plot spike
            errorbar(p+0.1,dat_not(p),...
                iqr_not(p,2) - dat_not(p),dat_not(p) - iqr_not(p,1),'ko',...
                'markersize',15,'linewidth',2);



        end


        %% Formatting
        if f == 2 || strcmp(met,'sd')
            xlabel('Patient');
        end
        xticklabels([])
        if f == 1
            legend('Spike','No spike','fontsize',20,'location','northwest')
        end
        if f == 1
            ylabel(sprintf('Pre-IED %s change',pretty_name));
        end
        if ~strcmp(met,'sd')
            title(sprintf(metrics.time.freq(f).name))
        end
        set(gca,'fontsize',20)

        %% Add p-value
        yl = get(gca,'ylim');
        xl = get(gca,'xlim');
        p_xloc = 0.5*(xl(1)+xl(2));
        p_yloc = yl(1) + 0.9*(yl(2)-yl(1));
        pp = pretty_p(pval,nfreq);
        text(p_xloc,p_yloc,pp,'horizontalalignment','center','fontsize',20);



    end
else
    figure

    if strcmp(met,'sd')
        nfreq = 1;
        set(gcf,'position',[1 100 400 270])
        [ha, pos] = tight_subplot(1, nfreq, [0.10 0.05], [0.12 0.11], [0.24 0.01]);
    else
        nfreq = length(metrics.time.freq);
        set(gcf,'position',[1 100 1100 270])
        [ha, pos] = tight_subplot(1, nfreq, [0.10 0.04], [0.12 0.11], [0.06 0.01]);
    end
    

for f = 1:nfreq
    axes(ha(f))
    dat_sp = metrics.time.freq(f).(met).auc.data(:,1);
    dat_not = metrics.time.freq(f).(met).auc.data(:,2);
    pval = metrics.time.freq(f).(met).auc.pval;
    np = length(dat_sp);
    x_sp = ones(np,1) + add_jitter(np,jitter_amount);
    x_not = 2*ones(np,1) + add_jitter(np,jitter_amount);
    
    plot(x_sp,dat_sp,'ro','markersize',10,'linewidth',2)
    hold on
    plot(x_not,dat_not,'ko','markersize',10,'linewidth',2)
    
    xticks([1 2])
    xticklabels({'IED','Not IED'})
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
    
end

print(gcf,[out_folder,sprintf('auc_%s',met)],'-dpng');

end

function out = add_jitter(n,amount)
    out = randn(n,1)*amount;

end