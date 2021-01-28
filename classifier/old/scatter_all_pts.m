function scatter_all_pts(stats,windows,met,which_pre_rise)
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

network_count = length(stats);
nfreq = length(stats(1).time(1).freq);

%% Initialize figure
figure

if skip_sd == 1
    
    set(gcf,'position',[1 100 900 280])
    
    [ha, pos] = tight_subplot(1, nfreq, [0.10 0.01], [0.12 0.11], [0.11 0.01]);
else
    set(gcf,'position',[1 100 900 550])
    [ha, pos] = tight_subplot(2, nfreq, [0.10 0.01], [0.10 0.06], [0.11 0.01]);
end

% Find the right window
for t = 1:length(stats(1).time)
    if ismember(stats(1).time(t).time_window,windows), break; end
end

%% First plot the SD
curr_sp = 0;
z_range = [0 0];
if skip_sd == 0

axes(ha(1))
set(ha(1),'Position',[pos{1}(1) pos{1}(2) pos{1}(3)*nfreq pos{1}(4)]);
dat_sp = stats(1).time(t).freq(1).sd.all_z_spike;
dat_not = stats(1).time(t).freq(1).sd.all_z_not;
times = stats(1).time(t).freq(1).sd.pt(1).times;


[z_range,prettyp] = do_scatter(dat_sp,dat_not,1,z_range);
xl = get(gca,'xlim');
yl = get(gca,'ylim');
text((xl(1)+xl(2))/2,yl(1)+0.9*(yl(2)-yl(1)),...
    sprintf('%s',prettyp),'fontsize',20,...
    'horizontalalignment','center')
ylabel('Normalized power');
title('Total power')

%% Delete next 2
for d = 2:nfreq
    delete(ha(d))
end
curr_sp = nfreq;
end

all_p = cell(nfreq,2);
pcount = 1;
%% Now plot the metric
for f = 1:nfreq
    sp = curr_sp+f;
    axes(ha(sp))
    
    dat_sp = stats(1).time(1).freq(f).(met).all_z_spike;
    dat_not = stats(1).time(1).freq(f).(met).all_z_not;
    
    [z_range,prettyp] = do_scatter(dat_sp,dat_not,nfreq,z_range);
    all_p{pcount,1} = prettyp;
    all_p{pcount,2} = sp;
    title(sprintf('%s',stats(1).time(1).freq(f).name))
    pcount = pcount + 1;
end

if skip_sd == 0
    curr_sp = nfreq;
else
    curr_sp = 0;
end

%% Change ylims and add p values
for sp = curr_sp+1:curr_sp+nfreq
    axes(ha(sp))
    
    % get ylim
    yl(1) = z_range(1) - 0.1*(z_range(2)-z_range(1));
    yl(2) = z_range(2) + 0.1*(z_range(2)-z_range(1));
    
    ylim(yl)
    xl = get(gca,'xlim');
    yl = get(gca,'ylim');
    if mod(sp,(nfreq)) ~= 1
        yticklabels([])
    end
    
    if mod(sp,nfreq) == 1
        ylabel(sprintf('Normalized%s',met_text))
    end
    
    for p = 1:length(all_p)
        if all_p{p,2} == sp
            text((xl(1)+xl(2))/2,yl(1)+0.9*(yl(2)-yl(1)),...
                sprintf('%s',all_p{p,1}),'fontsize',20,...
                'horizontalalignment','center')
        end
    end
end

%delete(ha(1))

%% Print
print(gcf,[out_folder,met,'_',num2str(which_pre_rise),'_scatter'],'-dpng');

end

function [z_range,prettyp] = do_scatter(dat_sp,dat_not,alpha,z_range)
% Remove columns with only one non nan
for j = 1:size(dat_sp,2)
    if sum(~isnan(dat_sp(:,j))) == 1
        dat_sp(:,j) = nan;
    end

    if sum(~isnan(dat_not(:,j))) == 1
        dat_not(:,j) = nan;
    end
end

% Take the last column of non nans
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

plot_sp = dat_sp(:,last_non_nan);   
plot_not = dat_not(:,last_non_nan);

hs = plot(randn(size(plot_sp,1),1)*0.05+ones(size(plot_sp,1),1)...
    ,plot_sp,'o','markersize',12);
set(hs, 'MarkerFaceColor', get(hs,'Color')); 
hold on
hn = plot(randn(size(plot_not,1),1)*0.05+2*ones(size(plot_not,1),1),...
    plot_not,'o','markersize',12);
set(hn, 'MarkerFaceColor', get(hn,'Color')); 
xlim([0.5 2.5])
xl = get(gca,'xlim');
xticks([1 2])
xticklabels({'Spike','No spike'})
set(gca,'fontsize',20)

%% Get range of values
min_data = min([plot_sp;plot_not]);
max_data = max([plot_sp;plot_not]);

if min_data < z_range(1)
    z_range(1) = min_data;
end

if max_data > z_range(2)
    z_range(2) = max_data;
end

end

