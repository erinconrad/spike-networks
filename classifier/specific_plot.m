function [z_range,prettyp] = specific_plot(dat_sp,dat_not,do_legend,z_range,adjust_z_range,times,alpha)

    %if isempty(dat_sp), z_range = z_range,prettyp = ''; return; end
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
    sp_pt = errorbar(times,dat_sp_mean,dat_sp_se,'linewidth',2);
    hold on
    not_pt = errorbar(times,dat_not_mean,dat_not_se,'linewidth',2);
    
    if do_legend == 1
        legend([sp_pt,not_pt],{'Spike','Spike-free'},'fontsize',20,...
        'location','northeastoutside');
    end
    
    %% Significance testing 
    % Test just the last point
    all_nan_columns = sum(isnan(dat_sp),1) == size(dat_sp,1);
    last_non_nan = find(all_nan_columns);
    last_non_nan(last_non_nan == 1) = [];
    if isempty(last_non_nan)
        last_non_nan = size(dat_sp,2);
    else
        last_non_nan = last_non_nan(1)-1;
    end
    [~,pval] = ttest(dat_sp(:,last_non_nan),dat_not(:,last_non_nan));
    prettyp = pretty_p(pval,alpha);


    maxyloc = [dat_sp_mean(~all_nan_columns)+dat_sp_se(~all_nan_columns),...
        dat_not_mean(~all_nan_columns)+dat_not_se(~all_nan_columns),...
    dat_sp_mean(~all_nan_columns)-dat_sp_se(~all_nan_columns),...
    dat_not_mean(~all_nan_columns)-dat_not_se(~all_nan_columns)];



    set(gca,'fontsize',20)
    
    if adjust_z_range == 1
        % adjust z_range
        if max(maxyloc) > z_range(2)
            z_range(2) = max(maxyloc);
        end

        if min(maxyloc) < z_range(1)
            z_range(1) = min(maxyloc);
        end
    end
    
    
end