function make_fig13(stats,windows,met,which_pre_rise)

if strcmp(met,'sd')
    met_text = 'power';
elseif strcmp(met,'ers')
    met_text = [newline,'frequency-specific power'];
elseif strcmp(met,'F')
    met_text = 'network difference';
else
    met_text = met;
end

locations = spike_network_files;
main_folder = locations.main_folder;
results_folder = [main_folder,'results/'];
out_folder = [results_folder,'plots/'];

if length(windows) > 1, error('Only one time window'); end

network_count = length(stats);
n_freq_abs = length(stats(1).time(1).freq);

%% Initialize figure
figure
set(gcf,'position',[1 100 1450 550])
[ha, pos] = tight_subplot(2, n_freq_abs, [0.15 0.01], [0.13 0.06], [0.07 0.01]);


z_range = [0 0];


% Find the right window
for t = 1:length(stats(1).time)
    if ismember(stats(1).time(t).time_window,windows), break; end
end


%% First plot the SD
axes(ha(1))
set(ha(1),'Position',[pos{1}(1) pos{1}(2) pos{1}(3)*4 pos{1}(4)]);
dat_sp = stats(1).time(t).freq(1).sd.all_z_spike;
dat_not = stats(1).time(t).freq(1).sd.all_z_not;
times = stats(1).time(t).freq(1).sd.pt(1).times;
[z_range,prettyp] = specific_plot(dat_sp,dat_not,1,z_range,0,times,1);
%xlabel('Time relative to spike peak (s)')
xl = get(gca,'xlim');
yl = get(gca,'ylim');
text((xl(1)+xl(2))/2,yl(1)+0.9*(yl(2)-yl(1)),...
    sprintf('%s',prettyp),'fontsize',20,...
    'horizontalalignment','center')
ylabel('Normalized power');

%% Delete next 6
for d = 2:n_freq_abs
    delete(ha(d))
end

all_p = cell((n_freq_abs+1)*2,2);
pcount = n_freq_abs+1;
%% Plot the ers
for n = 1:network_count

    %net_name = stats(n).name;
    tcount = 1;
    
    times = stats(n).time(t).freq(1).(met).pt(1).times;
    nfreq = length(stats(n).time(t).freq);
    
    for f = 1:nfreq
        
        % Get appropriate subplot
        sp = (n_freq_abs)*(tcount) + f;
        axes(ha(sp));
        
        dat_sp = stats(n).time(t).freq(f).(met).all_z_spike;
        dat_not = stats(n).time(t).freq(f).(met).all_z_not;
        
        % Plot them
        [z_range,prettyp] = specific_plot(dat_sp,dat_not,0,z_range,1,times,n_freq_abs);
        pcount = pcount + 1;
        all_p{pcount,2} = sp;
        all_p{pcount,1} = prettyp;
        
        if tcount == 1 
            title(sprintf('%s',...
                strrep(stats(n).time(t).freq(f).name,'_',' ')))
        end
        
        if (f == floor((n_freq_abs+1)/2) && n == 1) 
            xlabel('Time relative to spike peak (s)')
        end

        if n == network_count
            ylabel(sprintf('Normalized %s',met_text))    
        end

    end

end

for sp = n_freq_abs+1:length(ha)
    axes(ha(sp))
    %t = ceil(sp/(n_freq_abs+1));
    ylim(z_range)
    xl = get(gca,'xlim');
    yl = get(gca,'ylim');
    
    if mod(sp,(n_freq_abs)) ~= 1
        yticklabels([])
    end
    
    for p = 1:length(all_p)
        if all_p{p,2} == sp
            text((xl(1)+xl(2))/2,yl(1)+0.9*(yl(2)-yl(1)),...
                sprintf('%s',all_p{p,1}),'fontsize',20,...
                'horizontalalignment','center')
        end
    end
end

print(gcf,[out_folder,'Figure_',sprintf('%s_%d',met,which_pre_rise)],'-depsc');

end