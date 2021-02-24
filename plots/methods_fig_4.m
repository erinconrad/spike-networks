function methods_fig_4

%{
For network comparison
Panel 1: multiple spikes
Panel 2: 1 spike, time windows
Panel 2: Show adjacency matrices over time for single spike
Panel 3: vectorize matrices and compare vectors across times
%}

whichPt = 7;
which_times = 1:21;
windows = [-2:0.1:0];
time_window = 0.1;
which_spike = 6;
surround_time = 0.5;

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
seq_folder = [results_folder,'seq_data/'];

pt = load(pt_file); % will create a structure called "pt"
pt = pt.pt;
name = pt(whichPt).name;

spike = load([eeg_folder,sprintf('%s_eeg.mat',name)]);
spike = spike.spike;

seq = load([seq_folder,name,'_seq.mat']);
seq = seq.seq;

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

%% Single sequence
axes(ha(1))
s = which_spike;
values = spike(s).data;
fs = spike(s).fs;
data = pre_processing(values,1,0,1,fs);
curr_seq = seq(s).seq;
out_chs = curr_seq(:,1);
first = seq(s).first_ch;
biggest = seq(s).biggest_ch;
surround_idx = round(size(values,1)/2-fs*surround_time):round(size(values,1)/2+fs*surround_time);
offset = 0;
ch_bl = zeros(size(curr_seq,1),1);
for ich = 1:size(curr_seq,1)
    plot(linspace(-surround_time,surround_time,size(data(surround_idx,:),1)),data(surround_idx,out_chs(ich))+offset,'k')
    ch_bl(ich) = offset + median(data(surround_idx,out_chs(ich)));
    if ich<size(curr_seq,1)
        offset = offset - (max(data(surround_idx,out_chs(ich))) - min(data(surround_idx,out_chs(ich+1))));
    end
    hold on

end

yl = get(gca,'ylim');
xl = get(gca,'xlim');
arrow_x = [-0.1 -0.05];
arrow_y = [ch_bl(1) ch_bl(end)];

xlim([-surround_time-0.3,surround_time]);
text(-surround_time-0.3,ch_bl(out_chs==first),'Lead IED','fontsize',20)
text(-surround_time-0.3,ch_bl(out_chs==biggest),'Peak IED','fontsize',20)
myarrow(arrow_x,arrow_y);
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
plot(x_sp,dat_sp,'ro','markersize',10,'linewidth',2)
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
plot(x_sp,dat_sp,'ro','markersize',10,'linewidth',2)
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


print(gcf,[out_folder,'fig2'],'-depsc')

end

function myarrow(x,y)
ax = gca;
axpos = get(ax, 'Position');
X = get(gca,'XLim');
Y = get(gca,'YLim'); 
difX = X(2) - X(1);
difY = Y(2) - Y(1);
newx = (x-X(1))./difX;
newy = (y-Y(1))./difY;
annotation('arrow',[newx(1)*axpos(3)+axpos(1) newx(2)*axpos(3)+axpos(1)],...
    [newy(1)*axpos(4)+axpos(2) newy(2)*axpos(4)+axpos(2)],'linewidth',2)
end

function out = add_jitter(n,amount)
    out = randn(n,1)*amount;

end

