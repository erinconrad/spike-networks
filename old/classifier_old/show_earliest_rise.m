function show_earliest_rise(stats,windows,met,which_pre_rise)
skip_sd = 0;
if strcmp(met,'sd')
    met_text = 'power';
elseif strcmp(met,'ers')
    met_text = [newline,'frequency-specific power'];
elseif strcmp(met,'F')
    met_text = 'network difference';
elseif strcmp(met,'ge')
    met_text = ' global efficiency';
    skip_sd = 1;
elseif strcmp(met,'ns_big')
    met_text = [' spike electrode',newline,'node strength'];
    skip_sd = 1;
else
    met_text = met;
end

locations = spike_network_files;
main_folder = locations.main_folder;
results_folder = [main_folder,'results/'];
out_folder = [results_folder,'plots/'];

if length(windows) > 1, error('Only one time window'); end
nfreq = length(stats(1).time(1).freq);


%% Initialize figure
figure
set(gcf,'position',[1 100 1450 270])
if skip_sd == 1
    
    set(gcf,'position',[1 100 900 280])
    
    [ha, pos] = tight_subplot(1, nfreq, [0.10 0.01], [0.12 0.11], [0.11 0.01]);
else
    set(gcf,'position',[1 100 900 550])
    [ha, pos] = tight_subplot(2, nfreq, [0.10 0.01], [0.10 0.06], [0.11 0.01]);
end

z_range = [0 0];

% Find the right window
for t = 1:length(stats(1).time)
    if ismember(stats(1).time(t).time_window,windows), break; end
end
times = stats(1).time(t).freq(1).sd.pt(1).times;

%% First plot the SD
curr_sp = 0;
z_range = [0 0];
if skip_sd == 0

axes(ha(1))
set(ha(1),'Position',[pos{1}(1) pos{1}(2) pos{1}(3)*nfreq pos{1}(4)]);
dat_sp = stats(1).time(t).freq(1).sd.all_z_spike;
dat_not = stats(1).time(t).freq(1).sd.all_z_not;
times = stats(1).time(t).freq(1).sd.pt(1).times;

auc_sp = stats(1).time(t).freq(1).sd.auc_spike;
auc_not = stats(1).time(t).freq(1).sd.auc_not;


[z_range,prettyp,onep] = do_specific_plot(dat_sp,dat_not,0,z_range,1,times,1,auc_sp,auc_not);
xl = get(gca,'xlim');
yl = get(gca,'ylim');
text((xl(1)+xl(2))/2,yl(1)+0.9*(yl(2)-yl(1)),...
    sprintf('%s',onep),'fontsize',20,...
    'horizontalalignment','center')
ylabel('Normalized power');
title('Total power')

%% Delete next 2
for d = 2:nfreq
    delete(ha(d))
end
curr_sp = nfreq;
end

all_p = cell((nfreq),1);
all_one_p = zeros(nfreq,1);
%% Plot the ers 
for f = 1:nfreq
    sp = curr_sp+f;
    % Get appropriate subplot
    axes(ha(sp));

    dat_sp = stats(1).time(t).freq(f).(met).all_z_spike;
    dat_not = stats(1).time(t).freq(f).(met).all_z_not;
    
    auc_sp = stats(1).time(t).freq(f).(met).auc_spike;
    auc_not = stats(1).time(t).freq(f).(met).auc_not;

   % ts = stats.time.freq(f).(met).all_t_all_times;
    
    % Plot them
    [z_range,prettyp,one_p] = do_specific_plot(dat_sp,dat_not,0,z_range,1,times,nfreq,auc_sp,auc_not);
    all_p{f} = prettyp;
    all_one_p(f) = one_p;

    title(sprintf('%s',...
        strrep(stats(1).time(t).freq(f).name,'_',' ')))

    if (f == 2) 
        xlabel('Time relative to spike peak (s)')
    end

    if f == 1
        ylabel(sprintf('Normalized %s',met_text))    
    end

end



for f = 1:nfreq
    if ~skip_sd
        sp = nfreq+f;
    else
        sp = f;
    end
    axes(ha(sp))
    %t = ceil(sp/(n_freq_abs+1));
     % get ylim
    yl(1) = z_range(1) - 0.1*(z_range(2)-z_range(1));
    yl(2) = z_range(2) + 0.1*(z_range(2)-z_range(1));
    
    ylim(yl)
    
    xl = get(gca,'xlim');
    yl = get(gca,'ylim');
    
    if f ~= 1
        yticklabels([])
    end
    
    text((xl(1)+xl(2))/2,yl(1)+0.9*(yl(2)-yl(1)),...
                sprintf('%1.3f',all_one_p(f)),'fontsize',20,...
                'horizontalalignment','center')
    
    %{
    prettyp = all_p{sp};
    for t = 1:length(prettyp)
        text(times(t),yl(1)+0.9*(yl(2)-yl(1)),...
            sprintf('%s',prettyp{t}),'fontsize',30,...
                'horizontalalignment','center')
        %{
            text((xl(1)+xl(2))/2,yl(1)+0.9*(yl(2)-yl(1)),...
                sprintf('%s',all_p{p,1}),'fontsize',20,...
                'horizontalalignment','center')
        %}
        
    end
    %}
end



end



function [z_range,prettyp,one_p] = do_specific_plot(dat_sp,dat_not,do_legend,z_range,adjust_z_range,times,alpha,...
    auc_sp,auc_not)

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
    % Loop through all times
    
    prettyp = cell(size(dat_sp,2),1);
    for t = 1:size(dat_sp,2)
        [~,pval] = ttest(dat_sp(:,t),dat_not(:,t));
    	pp = get_asterisks(pval,alpha);
        prettyp{t} = pp;
    end
    %}
    
    %% Look at all times
    % For each patient, get a tstat from a paired ttest comparing the zs at
    % each time from spike and not spike
    %{
    ts = zeros(size(dat_sp,1),1);
    for p = 1:size(dat_sp,1)
        [~,~,~,st] = ttest(dat_sp(p,:),dat_not(p,:));
        ts(p) = st.tstat;
    end
    
    % one sample ttest of tstats
    
    [~,pval] = ttest(ts);
    one_p = pval;
    %}
    
    %% AUC
    % Take AUC for each pt
    [~,pv] = ttest(auc_sp,auc_not);
    one_p = pv;
    
   
    
    %% Do the spike metrics tend to be higher than the non-spike metrics?
    %{
    % Take the median across pts for each time
    median_dat_sp = median(dat_sp,1);
    median_dat_not = median(dat_not,1);
    
    % Paired ttest of sp and not sp for each time to see if they are higher
    % for the spike (remove first time as other times are relative to this)
    [~,pval] = ttest(median_dat_sp(2:end),median_dat_not(2:end));
    one_p = pval;
    %}
    
    %{
    new_ts = nan(size(ts,1),1);
    for p = 1:size(ts,1)
        [~,~,~,st] = ttest(ts(p,:));
        new_ts(p) = st.tstat;
    end
    %}
    
    % one sample ttest of tstats
    %{
    [~,pval] = ttest(new_ts);
    one_p = pval;
    %}
    
    maxyloc = [dat_sp_mean+dat_sp_se,...
        dat_not_mean+dat_not_se,...
    dat_sp_mean-dat_sp_se,...
    dat_not_mean-dat_not_se];
    %}



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