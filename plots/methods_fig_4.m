function methods_fig_4(metrics_sd)
metrics = metrics_sd;

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
set(gcf,'position',[1 300 900 500]);
[ha, pos] = tight_subplot(2, 2, [0.14 0.10], [0.07 0.02], [0.05 0.03]);


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
    hold on
    ch_bl(ich) = offset + median(data(surround_idx,out_chs(ich)));
    if curr_seq(ich) == first
        [~,max_idx] = min(data(surround_idx,out_chs(ich)));
    else
        [~,max_idx] = max(data(surround_idx,out_chs(ich)));
    end
    max_idx = max_idx + surround_idx(1)-1;
    max_time = (max_idx-1)/fs-3;
    if curr_seq(ich) == first
        fp = plot(max_time,data(max_idx,out_chs(ich))+offset,'o',...
            'markersize',10,'markerfacecolor',[0, 0.4470, 0.7410],'markeredgecolor','k');
    elseif curr_seq(ich) == biggest
        bp = plot(max_time,data(max_idx,out_chs(ich))+offset,'d',...
            'markersize',10,'markerfacecolor',[0.8500, 0.3250, 0.0980],'markeredgecolor','k');
    end
    if ich<size(curr_seq,1)
        offset = offset - 1.2*(max(data(surround_idx,out_chs(ich))) - min(data(surround_idx,out_chs(ich+1))));
    end
    

end

yl = get(gca,'ylim');
xl = get(gca,'xlim');
arrow_x = [0.08 0.14];
arrow_y = [ch_bl(1) ch_bl(end)];
lp=legend([fp,bp],{'Lead electrode','Peak electrode'},'fontsize',20,'location','southeast');
set(lp,'position',[0.3144 0.5375 0.1956 0.1090]);
%{
text(-surround_time-0.3,ch_bl(out_chs==first),'Lead IED','fontsize',20)
text(-surround_time-0.3,ch_bl(out_chs==biggest),'Peak IED','fontsize',20)
%}
myarrow(arrow_x,arrow_y);
ha(1).Visible = 'off';
set(findall(gca, 'type', 'text'), 'visible', 'on')

%% Compare SOZ for peak and lead IEDs
[percs,all] = soz_info;
axes(ha(2))
bar(percs*100);
%legend({'Lead electrode','Peak electrode'},'fontsize',20,'location','northeast')
xlabel('Patient')
ylabel('% in SOZ');
set(gca,'fontsize',20)


%% IEDs in SOZ compared to IEDs not in SOZ
axes(ha(3))
dat_soz = metrics.time.freq(1).sd_auto.auc.soz.data(:,1);
dat_not = metrics.time.freq(1).sd_auto.auc.soz.data(:,2);
pval = metrics.time.freq(1).sd_auto.auc.soz.pval;
np = length(dat_soz);
x_soz = ones(np,1) + add_jitter(np,0.05);
x_not = 2*ones(np,1) + add_jitter(np,0.05);

plot(x_soz,dat_soz,'ro','markersize',10,'linewidth',2)
hold on
plot(x_not,dat_not,'ko','markersize',10,'linewidth',2)

xticks([1 2])
xticklabels({'SOZ','Not SOZ'})

ylabel(sprintf('Relative power'))

xlim([0.5 2.5])
set(gca,'fontsize',20)

p_pretty = get_asterisks(pval,1);
yl = get(gca,'ylim');
set(gca,'ylim',[yl(1) yl(1) + 1.18*(yl(2)-yl(1))])
yl = get(gca,'ylim');
p_xloc = 1.5;
p_yloc = yl(1) + 0.95*(yl(2)-yl(1));
line_yloc =  yl(1) + 0.90*(yl(2)-yl(1));
text(p_xloc,p_yloc,p_pretty,'HorizontalAlignment','Center','fontsize',20)
plot([1 2],[line_yloc line_yloc],'k')


%% Lead IED electrode vs other electrodes in sequence
axes(ha(4))
dat_lead = metrics.time.freq(1).sd_auto.auc.first_v_other.data(:,1);
dat_other = metrics.time.freq(1).sd_auto.auc.first_v_other.data(:,2);
pval = metrics.time.freq(1).sd_auto.auc.first_v_other.pval;
np = length(dat_lead);
x_lead = ones(np,1) + add_jitter(np,0.05);
x_other = 2*ones(np,1) + add_jitter(np,0.05);

plot(x_lead,dat_lead,'ro','markersize',10,'linewidth',2)
hold on
plot(x_other,dat_other,'ko','markersize',10,'linewidth',2)

xticks([1 2])
xticklabels({'Lead electrode','Non-lead'})
ylabel(sprintf('Relative power'))
xlim([0.5 2.5])
set(gca,'fontsize',20)

p_pretty = get_asterisks(pval,1);
yl = get(gca,'ylim');
set(gca,'ylim',[yl(1) yl(1) + 1.18*(yl(2)-yl(1))])
yl = get(gca,'ylim');
p_xloc = 1.5;
p_yloc = yl(1) + 0.95*(yl(2)-yl(1));
line_yloc =  yl(1) + 0.90*(yl(2)-yl(1));
text(p_xloc,p_yloc,p_pretty,'HorizontalAlignment','Center','fontsize',20)
plot([1 2],[line_yloc line_yloc],'k')

%% annotations
annotation('textbox',[0.01 0.92 0.1 0.1],'String','A','linestyle','none','fontsize',30);
annotation('textbox',[0.48 0.92 0.1 0.1],'String','B','linestyle','none','fontsize',30);
annotation('textbox',[0.01 0.40 0.1 0.1],'String','C','linestyle','none','fontsize',30);
annotation('textbox',[0.48 0.40 0.1 0.1],'String','D','linestyle','none','fontsize',30);

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
    [newy(1)*axpos(4)+axpos(2) newy(2)*axpos(4)+axpos(2)],'linewidth',3)
end

function out = add_jitter(n,amount)
    out = randn(n,1)*amount;

end

