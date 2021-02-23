function methods_fig_2(metrics_sd,metrics_ers,earliest_rise,orig_pt_rise)

%{
For network comparison
Panel 1: multiple spikes
Panel 2: 1 spike, time windows
Panel 2: Show adjacency matrices over time for single spike
Panel 3: vectorize matrices and compare vectors across times
%}

whichPt = 7;
surround = 2;
full = 3;
time_window = 0.1;
whichSpikes = [3 4 7];
which_times = 1:12;
windows = [-2:0.1:0];

%% Get file locations, load spike times and pt structure
locations = spike_network_files;
main_folder = locations.main_folder;
eeg_folder = [main_folder,'results/eeg_data/'];
results_folder = [main_folder,'results/'];
data_folder = [main_folder,'data/'];
script_folder = locations.script_folder;
addpath(genpath(script_folder));
pt_file = [data_folder,'spike_structures/pt.mat'];
out_folder = [results_folder,'plots/'];

pt = load(pt_file); % will create a structure called "pt"
pt = pt.pt;
name = pt(whichPt).name;

spike = load([eeg_folder,sprintf('%s_eeg.mat',name)]);
spike = spike.spike;



figure
set(gcf,'position',[1 300 900 750]);
[ha, pos] = tight_subplot(3, 3, [0.08 0], [0.05 0.05], [0.06 0.02]);
xgap = 0.005;
delete(ha(3));
delete(ha(6));
set(ha(1),'position',[pos{1}(1) pos{1}(2) pos{1}(3)*1.5-xgap pos{1}(4)])
set(ha(2),'position',[pos{1}(1)+pos{1}(3)*1.5+xgap pos{2}(2) pos{1}(3)*1.5-xgap pos{2}(4)])

set(ha(4),'position',[pos{4}(1) pos{4}(2) pos{4}(3)*1.5-xgap pos{4}(4)])
set(ha(5),'position',[pos{4}(1)+pos{4}(3)*1.5+xgap pos{5}(2) pos{4}(3)*1.5-xgap pos{5}(4)])

set(ha(7),'position',[pos{7}(1) pos{7}(2) pos{7}(3)-xgap pos{7}(4)])
set(ha(8),'position',[pos{8}(1)+xgap pos{8}(2) pos{8}(3)-xgap pos{8}(4)])
set(ha(9),'position',[pos{9}(1)+2*xgap pos{9}(2) pos{9}(3)-xgap pos{9}(4)])


%% Example spikes
axes(ha(1))
offset = 0;
time_offset = 0;
count = 0;
for s = whichSpikes
    count = count + 1;
    ch = spike(s).biggest_dev;
    values = spike(s).data;
    fs = spike(s).fs;
    values = pre_processing(values,1,0,1,fs);
    sindex = (full-surround)*fs;
    sindex = [round(sindex),round(size(values,1)-sindex)];
    eeg = values(sindex(1):sindex(2),ch);
    times = linspace(-surround+time_offset,surround+time_offset,length(eeg));
    plot(times,eeg+offset,'k','linewidth',2)
    if count < length(whichSpikes)
        offset = offset + ...
            0.7*abs(max(spike(whichSpikes(count+1)).data(:,spike(whichSpikes(count+1)).biggest_dev) ...
            - min(eeg)));
        time_offset = time_offset + 1;
    end
    hold on
    
    % Plot dotte
end
%xlabel('Time (s) relative to spike peak')
xticklabels([])
yticklabels([])
%{
text(0.5,1,'Manually identify spikes','Units','Normalized',...
    'HorizontalAlignment','Center','Fontsize',20)
%}
%title('Manually identify spikes')
set(gca,'fontsize',20)
ha(1).Visible = 'off';
set(findall(gca, 'type', 'text'), 'visible', 'on')

%% Single spike with time windows
axes(ha(2))
s = whichSpikes(2);
ch = spike(s).biggest_dev;
values = spike(s).data;
values = pre_processing(values,1,0,1,fs);
values = values(:,ch);
if whichPt == 7
    new_pt = 2;
else
    error('need to assign correct pt');
end
rise = min(orig_pt_rise(new_pt).both(s,:));

% Get index windows
peak = round(length(values)/2);
n_windows = length(windows);
index_windows = zeros(n_windows,2);

for i = 1:n_windows
    index_windows(i,1) = peak + round(windows(i)*fs);
    index_windows(i,2) = peak + round(windows(i)*fs) + round((windows(2)-windows(1))*fs);
end
values = values(index_windows(1,1):index_windows(end,2));
times = linspace(windows(1),windows(end)+0.1,length(values));
plot(times,values,'k','linewidth',2);
line_times = -2:0.1:0.1;
hold on
for i = 1:length(line_times)
    if line_times(i) + 0.1 >= rise
        pe = plot([line_times(i) line_times(i)],get(ha(2),'ylim'),'r:','linewidth',2);
    else
        plot([line_times(i) line_times(i)],get(ha(2),'ylim'),'k--','linewidth',2);
    end
    if i < length(line_times)
    end
end
%xlabel('Time (s) relative to spike peak')
legend(pe,'Excluded windows','fontsize',20,'location','southwest')
yticklabels([])
%{
text(0.5,1,'Divide EEG signal into time windows','Units','Normalized',...
    'HorizontalAlignment','Center')
%}
%title('Divide EEG signal into time windows')
set(gca,'fontsize',20)
xlim([-2 0.1])
ha(2).Visible = 'off';
set(findall(gca, 'type', 'text'), 'visible', 'on')

%% AUC
metrics = metrics_sd;
axes(ha(4))
met = 'sd_auto';
f = 1;
dat_sp = metrics.time.freq(f).(met).auc.data(:,1);
dat_not = metrics.time.freq(f).(met).auc.data(:,2);
pval = metrics.time.freq(f).(met).auc.pval;
np = length(dat_sp);
x_sp = ones(np,1) + add_jitter(np,0.05);
x_not = 2*ones(np,1) + add_jitter(np,0.05);
plot(x_sp,dat_sp,'ro','markersize',10,'linewidth',2)
hold on
plot(x_not,dat_not,'ko','markersize',10,'linewidth',2)
xticks([1 2])
xticklabels({'IED','Not IED'})
ylabel({'Relative power'});
xlim([0.5 2.5])
set(gca,'fontsize',20)
p_pretty = get_asterisks(pval,1);
yl = get(gca,'ylim');
p_xloc = 1.5;
p_yloc = yl(1) + 1.1*(yl(2)-yl(1));
text(p_xloc,p_yloc,p_pretty,'HorizontalAlignment','Center','fontsize',30)
line_yloc =  yl(1) + 1.07*(yl(2)-yl(1));
plot([1 2],[line_yloc line_yloc],'k')
set(gca,'ylim',[yl(1) yl(1) + 1.18*(yl(2)-yl(1))])
yl = get(gca,'ylim');

%% Time change
axes(ha(5))
data = metrics.time.freq(f).(met).auc.short.data;
pvals = metrics.time.freq(f).(met).auc.short.ps;
ntimes = size(data,1);
times = metrics.time.freq(f).(met).auc.times;

min_rise = min(earliest_rise,[],2); % earliest reviewer notation for each spike
mean_rise_spikes = mean(min_rise); % mean across spikes
before_rise = zeros(ntimes,1);
last_before_rise = 1;
% Round to lowest 0.1 s
for tt = 1:ntimes
    if mean_rise_spikes > times(tt)
        before_rise(tt) = 1;
        last_before_rise = tt;
    end
end


auc_spike = squeeze(data(:,:,1));
mean_auc_spike = nanmean(auc_spike,2);
std_spike = nanstd(auc_spike,0,2);

auc_not = squeeze(data(:,:,2));
mean_auc_not = nanmean(auc_not,2);
std_not = nanstd(auc_not,0,2);

errorbar(times(1:last_before_rise-1),...
    mean_auc_spike(1:last_before_rise-1)...
    ,std_spike(1:last_before_rise-1),'ro','markersize',15,...
    'linewidth',2)

hold on

errorbar(times(1:last_before_rise-1),...
    mean_auc_not(1:last_before_rise-1)...
    ,std_not(1:last_before_rise-1),'ko','markersize',15,...
    'linewidth',2)

%xlabel('Time (s)');
set(gca,'fontsize',20)
labels = get(gca,'xticklabels');
new_labels = cell(length(labels),1);
for i = 1:length(labels)
    new_labels{i} = [labels{i},' s'];
end
xticklabels(new_labels)

sub_alpha = zeros(ntimes,1);
for tt = 1:last_before_rise-1
    num_left = last_before_rise-1-tt+1;
    num_sub_alpha = sum(pvals(tt:last_before_rise-1) < 0.05);
    if num_sub_alpha == num_left
        sub_alpha(tt) = 1;
    end
end

ylim(yl);
%xlim([-2 0.1])
yl = get(gca,'ylim');
p_yloc = yl(1) + 0.9*(yl(2)-yl(1));
for tt = 1:ntimes
    if sub_alpha(tt) == 1
        text(times(tt),p_yloc,'*','horizontalalignment','center','fontsize',30);
    end
end
yticklabels([])

endh = plot([mean_rise_spikes mean_rise_spikes],get(gca,'ylim'),'k--','linewidth',2);

legend('IED','Not IED','Mean visual rise time','fontsize',20,'location','northwest')

%% ERS
ps = cell(3,1);
metrics = metrics_ers;
y_range = [0 0];
for f = 1:3
axes(ha(6+f))
met = 'ers_auto';
dat_sp = metrics.time.freq(f).(met).auc.data(:,1);
dat_not = metrics.time.freq(f).(met).auc.data(:,2);
pval = metrics.time.freq(f).(met).auc.pval;
np = length(dat_sp);
x_sp = ones(np,1) + add_jitter(np,0.05);
x_not = 2*ones(np,1) + add_jitter(np,0.05);
plot(x_sp,dat_sp,'ro','markersize',10,'linewidth',2)
hold on
plot(x_not,dat_not,'ko','markersize',10,'linewidth',2)
xticks([1 2])
xticklabels({'IED','Not IED'})
if f == 1
    ylabel({'Relative power'});
end
title(metrics.time.freq(f).name)
xlim([0.5 2.5])
set(gca,'fontsize',20)
p_pretty = get_asterisks(pval,1);
ps{f} = p_pretty;
yl = get(gca,'ylim');
set(gca,'ylim',[yl(1) yl(1) + 1.18*(yl(2)-yl(1))])
yl = get(gca,'ylim');
if yl(1) < y_range(1)
    y_range(1) = yl(1);
end

if yl(2) > y_range(2)
    y_range(2) = yl(2);
end

if f ~= 1
    yticklabels([])
end

end

for f =1:3
    axes(ha(6+f))
    ylim([y_range(1) y_range(2)])
    yl = get(gca,'ylim');
    p_xloc = 1.5;
    p_yloc = yl(1) + 0.92*(yl(2)-yl(1));
    text(p_xloc,p_yloc,p_pretty,'HorizontalAlignment','Center','fontsize',30)
    line_yloc =  yl(1) + 0.90*(yl(2)-yl(1));
    plot([1 2],[line_yloc line_yloc],'k')
    
end

%% annotations
annotation('textbox',[0.03 0.90 0.1 0.1],'String','A','linestyle','none','fontsize',30);
annotation('textbox',[0.52 0.90 0.1 0.1],'String','B','linestyle','none','fontsize',30);
annotation('textbox',[0.05 0.57 0.1 0.1],'String','C','linestyle','none','fontsize',30);
annotation('textbox',[0.52 0.57 0.1 0.1],'String','D','linestyle','none','fontsize',30);
annotation('textbox',[0.05 0.24 0.1 0.1],'String','E','linestyle','none','fontsize',30);
annotation('textbox',[0.37 0.24 0.1 0.1],'String','F','linestyle','none','fontsize',30);
annotation('textbox',[0.68 0.24 0.1 0.1],'String','G','linestyle','none','fontsize',30);

annotation('line',[0.41 0.52],[0.84 0.87])
annotation('line',[0.41 0.52],[0.81 0.75])

print(gcf,[out_folder,'fig1'],'-depsc')

end

function out = add_jitter(n,amount)
    out = randn(n,1)*amount;

end
