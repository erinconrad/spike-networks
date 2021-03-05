function methods_fig_3(metrics_ns_big,metrics_ns_avg,earliest_rise,orig_pt_rise)

%{
For network comparison
Panel 1: multiple spikes
Panel 2: 1 spike, time windows
Panel 2: Show adjacency matrices over time for single spike
Panel 3: vectorize matrices and compare vectors across times
%}

whichPt = 7;
duration = 6;
whichSpikes = [3 4 7];
which_times = 1:21;
windows = [-2:0.1:0];
time_window = 0.1;

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

%
adj_folder = [results_folder,'adj_mat/manual/adj_coherence/',sprintf('%1.1f/',time_window)];
meta = load([adj_folder,name,'_adj.mat']);
meta = meta.meta;
%}

figure
set(gcf,'position',[1 300 900 750]);
[ha, pos] = tight_subplot(3, 3, [0.08 0], [0.05 0.05], [0.12 0.03]);
xgap = 0.005;
delete(ha(3));
set(ha(1),'position',[pos{1}(1) pos{1}(2) pos{1}(3)*1.5-xgap pos{1}(4)])
set(ha(2),'position',[pos{1}(1)+pos{1}(3)*1.5+xgap pos{2}(2) pos{1}(3)*1.5-xgap pos{2}(4)])

set(ha(4),'position',[pos{4}(1) pos{4}(2) pos{4}(3)-xgap pos{4}(4)])
set(ha(5),'position',[pos{5}(1)+xgap pos{5}(2) pos{5}(3)-xgap pos{5}(4)])
set(ha(6),'position',[pos{6}(1)+2*xgap pos{6}(2) pos{6}(3)-xgap pos{6}(4)])


set(ha(7),'position',[pos{7}(1) pos{7}(2) pos{7}(3)-xgap pos{7}(4)])
set(ha(8),'position',[pos{8}(1)+xgap pos{8}(2) pos{8}(3)-xgap pos{8}(4)])
set(ha(9),'position',[pos{9}(1)+2*xgap pos{9}(2) pos{9}(3)-xgap pos{9}(4)])

% Single spike with time windows
axes(ha(1))
s = whichSpikes(2);
ch = spike(s).biggest_dev;
values = spike(s).data;
fs = spike(s).fs;
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
        pe = plot([line_times(i) line_times(i)],get(ha(1),'ylim'),'r:','linewidth',2);
    else
        plot([line_times(i) line_times(i)],get(ha(1),'ylim'),'k--','linewidth',2);
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
ha(1).Visible = 'off';
set(findall(gca, 'type', 'text'), 'visible', 'on')

% Adjacency matrices at multiple times
axes(ha(2))
s = whichSpikes(3);
adj_all_t = meta.spike(s).adj.adj;
y_offset = 0;
x_offset = 0;
for t = 1:length(which_times)
    imagesc([x_offset x_offset+size(adj_all_t,2)],...
        [y_offset y_offset+size(adj_all_t,3)],...
        squeeze(adj_all_t(which_times(t),:,:)));
    hold on
    plot([x_offset x_offset],...
        [y_offset y_offset+size(adj_all_t,3)],...
        'k','linewidth',2);
    plot([x_offset+size(adj_all_t,2) x_offset+size(adj_all_t,2)],...
        [y_offset y_offset+size(adj_all_t,3)],...
        'k','linewidth',2);
    plot([x_offset x_offset+size(adj_all_t,2)],...
        [y_offset y_offset],...
        'k','linewidth',2);
    plot([x_offset x_offset+size(adj_all_t,2)],...
        [y_offset+size(adj_all_t,3) y_offset+size(adj_all_t,3)],...
        'k','linewidth',2);
    if t == 2
        arrow_start_x = x_offset;
        arrow_start_y = y_offset+size(adj_all_t,3)+20;
    elseif t == 20
        arrow_end_x = x_offset;
        arrow_end_y = y_offset+size(adj_all_t,3)+20;
        
    elseif t == 21
        x_ch_arrow_l = [x_offset x_offset+size(adj_all_t,2)];
        y_ch_arrow_l = [y_offset+size(adj_all_t,3) y_offset+size(adj_all_t,3)];

        
        x_ch_arrow_r = [x_offset+size(adj_all_t,2) x_offset+size(adj_all_t,2)];
        y_ch_arrow_r = [y_offset+size(adj_all_t,3) y_offset];
    end
    if t < length(which_times)
        y_offset = y_offset + size(adj_all_t,2)*0.2;
        x_offset = x_offset + size(adj_all_t,2)*0.2;
    end
    
end
xlim([0 y_offset+size(adj_all_t,2)])
ylim([0 y_offset+size(adj_all_t,2)])
y = [arrow_start_y arrow_end_y];
yl = get(gca,'ylim');
y = yl(2) - y;
myarrow([arrow_start_x arrow_end_x],y);

axpos = get(gca, 'Position');

text(mean(x_ch_arrow_l),mean(y_ch_arrow_l)+38,'Electrode #',...
    'fontsize',17,'HorizontalAlignment','Center')
ang2 = atan2((y_ch_arrow_r(2) - y_ch_arrow_r(1))*axpos(4), (x_ch_arrow_r(2) - x_ch_arrow_r(1))*axpos(3)) * 180 / pi;
text(mean(x_ch_arrow_r)+21,mean(y_ch_arrow_r),'Electrode #',...
    'fontsize',17,'HorizontalAlignment','Center','rotation',-ang2)

y_ch_arrow_l = yl(2)-y_ch_arrow_l;
myarrow(x_ch_arrow_l,y_ch_arrow_l-18);
y_ch_arrow_r = yl(2)-y_ch_arrow_r;
myarrow(x_ch_arrow_r+9,y_ch_arrow_r);

%X = get(gca,'XLim');
%Y = get(gca,'YLim'); 
%difX = X(2) - X(1);
%difY = Y(2) - Y(1);


ang = atan2((arrow_end_y - arrow_start_y)*axpos(4), (arrow_end_x - arrow_start_x)*axpos(3)) * 180 / pi;

text(mean([arrow_start_x,arrow_end_x]), mean([arrow_start_y arrow_end_y])+50,...
    'Time','fontsize',25,'rotation',-ang,'HorizontalAlignment','Center')



%title('Calculate network for each time window')
set(gca,'fontsize',20)
ha(2).Visible = 'Off';
set(findall(gca, 'type', 'text'), 'visible', 'on')

%% NS big
ps = cell(3,1);
metrics = metrics_ns_big;
y_range = [0 0];
for f = 1:3
axes(ha(3+f))
met = 'ns_auto';
dat_sp = metrics.time.freq(f).(met).auc.data(:,1);
dat_not = metrics.time.freq(f).(met).auc.data(:,2);
pval = metrics.time.freq(f).(met).auc.pval;
np = length(dat_sp);
x_sp = ones(np,1) + add_jitter(np,0.05);
x_not = 2*ones(np,1) + add_jitter(np,0.05);
plot(x_sp,dat_sp,'o','markersize',10,'linewidth',2)
hold on
plot(x_not,dat_not,'ko','markersize',10,'linewidth',2)
xticks([1 2])
xticklabels({'IED','Not IED'})
if f == 1
    ylabel({'Relative node strength','(peak IED electrode)'});
end
title(metrics.time.freq(f).name)
xlim([0.5 2.5])
set(gca,'fontsize',20)
p_pretty = get_asterisks(pval,3);
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
    axes(ha(3+f))
    ylim([y_range(1) y_range(2)])
    yl = get(gca,'ylim');
    p_xloc = 1.5;
    
    if strcmp(ps{f},'ns')
        p_yloc = yl(1) + 0.95*(yl(2)-yl(1));
        text(p_xloc,p_yloc,ps{f},'HorizontalAlignment','Center','fontsize',20)
    else
        p_yloc = yl(1) + 0.92*(yl(2)-yl(1));
        text(p_xloc,p_yloc,ps{f},'HorizontalAlignment','Center','fontsize',30)
    end
    line_yloc =  yl(1) + 0.90*(yl(2)-yl(1));
    plot([1 2],[line_yloc line_yloc],'k')
    
end

%% NS avg
ps = cell(3,1);
metrics = metrics_ns_avg;
y_range = [0 0];
for f = 1:3
axes(ha(6+f))
met = 'ns_avg';
dat_sp = metrics.time.freq(f).(met).auc.data(:,1);
dat_not = metrics.time.freq(f).(met).auc.data(:,2);
pval = metrics.time.freq(f).(met).auc.pval;
np = length(dat_sp);
x_sp = ones(np,1) + add_jitter(np,0.05);
x_not = 2*ones(np,1) + add_jitter(np,0.05);
plot(x_sp,dat_sp,'o','markersize',10,'linewidth',2)
hold on
plot(x_not,dat_not,'ko','markersize',10,'linewidth',2)
xticks([1 2])
xticklabels({'IED','Not IED'})
if f == 1
    ylabel({'Relative node strength','(electrode average)'});
end
%title(metrics.time.freq(f).name)
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
    
    if strcmp(ps{f},'ns')
        p_yloc = yl(1) + 0.95*(yl(2)-yl(1));
        text(p_xloc,p_yloc,ps{f},'HorizontalAlignment','Center','fontsize',20)
    else
        p_yloc = yl(1) + 0.92*(yl(2)-yl(1));
        text(p_xloc,p_yloc,ps{f},'HorizontalAlignment','Center','fontsize',30)
    end
    line_yloc =  yl(1) + 0.90*(yl(2)-yl(1));
    plot([1 2],[line_yloc line_yloc],'k')
    
end

%% annotations
annotation('textbox',[0.11 0.90 0.1 0.1],'String','A','linestyle','none','fontsize',30);
annotation('textbox',[0.54 0.90 0.1 0.1],'String','B','linestyle','none','fontsize',30);
annotation('textbox',[0.11 0.57 0.1 0.1],'String','C','linestyle','none','fontsize',30);
annotation('textbox',[0.40 0.57 0.1 0.1],'String','D','linestyle','none','fontsize',30);
annotation('textbox',[0.69 0.57 0.1 0.1],'String','E','linestyle','none','fontsize',30);
annotation('textbox',[0.11 0.24 0.1 0.1],'String','F','linestyle','none','fontsize',30);
annotation('textbox',[0.40 0.24 0.1 0.1],'String','G','linestyle','none','fontsize',30);
annotation('textbox',[0.69 0.24 0.1 0.1],'String','H','linestyle','none','fontsize',30);



print(gcf,[out_folder,'fig2'],'-depsc')

end

function myarrow(x,y)
ax = gca;
axpos = get(ax, 'Position');
X = get(gca,'XLim');
Y = get(gca,'YLim'); 
difX = X(2) - X(1);
difY = Y(2) - Y(1);
newx = x./difX;
newy = y./difY;
annotation('arrow',[newx(1)*axpos(3)+axpos(1) newx(2)*axpos(3)+axpos(1)],[newy(1)*axpos(4)+axpos(2) newy(2)*axpos(4)+axpos(2)])
end

function out = add_jitter(n,amount)
    out = randn(n,1)*amount;

end

