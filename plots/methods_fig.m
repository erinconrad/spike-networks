function methods_fig

%{
For network comparison
Panel 1: multiple spikes
Panel 2: 1 spike, time windows
Panel 2: Show adjacency matrices over time for single spike
Panel 3: vectorize matrices and compare vectors across times
%}

whichPt = 6;
duration = 6;
time_window = 0.5;
whichSpikes = [10 11 14];
which_times = 1:12;

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

adj_folder = [results_folder,'adj_mat/manual/adj_simple/',sprintf('%1.1f/',time_window)];
meta = load([adj_folder,name,'_adj.mat']);
meta = meta.meta;

figure
set(gcf,'position',[1 300 900 800]);
[ha, pos] = tight_subplot(2, 2, [0.15 0.01], [0.05 0.05], [0.03 0]);

% Example spikes
axes(ha(1))
offset = 0;
time_offset = 0;
count = 0;
for s = whichSpikes
    count = count + 1;
    ch = spike(s).biggest_dev;
    values = spike(s).data(:,ch);
    times = linspace(-3+time_offset,3+time_offset,length(values));
    plot(times,values+offset,'k','linewidth',2)
    if count < length(whichSpikes)
        offset = offset + ...
            abs(max(spike(whichSpikes(count+1)).data(:,spike(whichSpikes(count+1)).biggest_dev) ...
            - min(values)));
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
title('Manually identify spikes')
set(gca,'fontsize',20)
ha(1).Visible = 'off';
set(findall(gca, 'type', 'text'), 'visible', 'on')

% Single spike with time windows
axes(ha(2))
s = whichSpikes(3);
ch = spike(s).biggest_dev;
values = spike(s).data(:,ch);
times = linspace(-3,3,length(values));
plot(times,values,'k','linewidth',2);
line_times = -duration/2:time_window:duration/2;
hold on
for i = 1:length(line_times)
    plot([line_times(i) line_times(i)],get(ha(2),'ylim'),'b--','linewidth',2);
    if i < length(line_times)
    end
end
%xlabel('Time (s) relative to spike peak')
yticklabels([])
%{
text(0.5,1,'Divide EEG signal into time windows','Units','Normalized',...
    'HorizontalAlignment','Center')
%}
title('Divide EEG signal into time windows')
set(gca,'fontsize',20)
ha(2).Visible = 'off';
set(findall(gca, 'type', 'text'), 'visible', 'on')

% Adjacency matrices at multiple times
axes(ha(4))
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
    elseif t == 11
        arrow_end_x = x_offset;
        arrow_end_y = y_offset+size(adj_all_t,3)+20;
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

%X = get(gca,'XLim');
%Y = get(gca,'YLim'); 
%difX = X(2) - X(1);
%difY = Y(2) - Y(1);
axpos = get(gca, 'Position');

ang = atan2((arrow_end_y - arrow_start_y)*axpos(4), (arrow_end_x - arrow_start_x)*axpos(3)) * 180 / pi;

text(mean([arrow_start_x,arrow_end_x]), mean([arrow_start_y arrow_end_y])+30,...
    'Time','fontsize',25,'rotation',-ang,'HorizontalAlignment','Center')
title('Calculate network for each time window')
set(gca,'fontsize',20)
ha(4).Visible = 'Off';
set(findall(gca, 'type', 'text'), 'visible', 'on')

% 2 times, all spikes, vectorized
n_spikes = length(spike);
nch = size(adj_all_t,2);
n_times = size(adj_all_t,1);
all_spikes = nan(n_spikes,nch*(nch-1)/2,n_times);

for s = 1:n_spikes
    for t = 1:size(meta.spike(s).adj.adj,1)
        adj_t = squeeze(meta.spike(s).adj.adj(t,:,:));
        vec_adj = flatten_or_expand_adj(adj_t);
        all_spikes(s,:,t) = vec_adj;
    end 
end

axes(ha(3))
initial_offset_y = size(all_spikes,1)*1.2;
imagesc([0 size(all_spikes,1)*1.5],...
    [initial_offset_y initial_offset_y+size(all_spikes,1)*1.5],all_spikes(:,:,1))
hold on
plot([0 0],[initial_offset_y initial_offset_y+size(all_spikes,1)*1.5],...
    'k','linewidth',2);
plot([0 size(all_spikes,1)*1.5],...
    [initial_offset_y+size(all_spikes,1)*1.5 initial_offset_y+size(all_spikes,1)*1.5],...
    'k','linewidth',2);
plot([size(all_spikes,1)*1.5 size(all_spikes,1)*1.5],...
    [initial_offset_y initial_offset_y+size(all_spikes,1)*1.5],...
    'k','linewidth',2);
plot([0 size(all_spikes,1)*1.5],[initial_offset_y initial_offset_y],...
    'k','linewidth',2);
text(mean([0 size(all_spikes,1)*1.5]),...
    initial_offset_y+size(all_spikes,1)*1.5+10,...
    'Time 1','FontSize',20,'HorizontalAlignment','Center');


text(mean([0 size(all_spikes,1)*1.5]),...
    initial_offset_y-15,...
    'Spike number','FontSize',20,'HorizontalAlignment','Center');
text(-15,...
    mean([initial_offset_y,initial_offset_y+size(all_spikes,1)*1.5]),...
    'Network element','FontSize',20,'HorizontalAlignment','Center',...
    'rotation',90);
x_offset = size(all_spikes,1)*2;
y_offset = 0;
for t = 2:size(all_spikes,3)
    imagesc([x_offset x_offset+size(all_spikes,1)],...
        [y_offset y_offset],all_spikes(:,:,t))
    
    plot([x_offset x_offset],...
        [y_offset y_offset+size(all_spikes,1)],...
        'k','linewidth',2);
    plot([x_offset x_offset+size(all_spikes,1)],...
        [y_offset+size(all_spikes,1) y_offset+size(all_spikes,1)],...
        'k','linewidth',2);
    plot([x_offset+size(all_spikes,1) x_offset+size(all_spikes,1)],...
        [y_offset y_offset+size(all_spikes,1)],...
        'k','linewidth',2);
    plot([x_offset x_offset+size(all_spikes,1)],...
        [y_offset y_offset],...
        'k','linewidth',2);
    
    
    x_offset = x_offset + size(all_spikes,1)*0.2;
    y_offset = y_offset + size(all_spikes,1)*0.2;
end
text(x_offset+ size(all_spikes,1)*.3,y_offset+ size(all_spikes,1),...
    'Other times','FontSize',20,'HorizontalAlignment','Center');
xlim([-10 x_offset+size(all_spikes,1)])
ylim([0 y_offset+size(all_spikes,1)])
title(sprintf('Compare networks for all spikes\nacross time windows'))
set(gca,'fontsize',20)


y = [initial_offset_y,initial_offset_y];
yl = get(gca,'ylim');
y = yl(2) - y;
myarrow([0 size(all_spikes,1)*1.5]+10,y+5)

y = [initial_offset_y,initial_offset_y+size(all_spikes,1)*1.5];
y = yl(2)-y;
myarrow([5 5],y);
xl = get(gca,'xlim');

text(mean(xl)-30,mean(yl),'vs.','fontsize',20);

ha(3).Visible = 'Off';    
set(findall(gca, 'type', 'text'), 'visible', 'on')

%% Inter-plot arrows, letters
annotation('arrow',[0.42 0.48],[0.78 0.78],'linewidth',2)
annotation('arrow',[0.7 0.7],[0.56 0.48],'linewidth',2)
annotation('arrow',[0.48 0.42],[0.3 0.3],'linewidth',2)
annotation('textbox',[0.04 0.9 0.1 0.1],'String','A',...
    'linestyle','none','Fontsize',30)
annotation('textbox',[0.47 0.9 0.1 0.1],'String','B',...
    'linestyle','none','Fontsize',30)
annotation('textbox',[0.47 0.43 0.1 0.1],'String','C',...
    'linestyle','none','Fontsize',30)
annotation('textbox',[0.04 0.43 0.1 0.1],'String','D',...
    'linestyle','none','Fontsize',30)
    

print([out_folder,'methods_fig'],'-depsc');

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

