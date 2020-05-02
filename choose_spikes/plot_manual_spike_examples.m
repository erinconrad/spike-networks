function plot_manual_spike_examples(avg)

%% Get locations
locations = spike_network_files;
script_folder = locations.script_folder;
main_folder = locations.main_folder;
pwname = locations.pwfile;
addpath(genpath(script_folder));

if isempty(locations.ieeg_folder) == 0
    addpath(genpath(locations.ieeg_folder));
end

pt_folder = [main_folder,'data/spike_structures/'];
pt = load([pt_folder,'pt.mat']);
pt = pt.pt;



eeg_folder = [main_folder,'results/eeg_data/'];
plot_folder = [eeg_folder,'sample_spikes/'];
if exist(plot_folder,'dir') == 0
    mkdir(plot_folder)
end

listing = dir([eeg_folder,'*.mat']);




for i = 1:length(listing)
    
    split_name = strsplit(listing(i).name,'_');
    ptname = split_name{1};
    
    eeg = load([eeg_folder,listing(i).name]);
    spike = eeg.spike;
    
    all_dev = [];
    
    if avg == 0
    
        for s = 1:length(spike)

            median_signal = median(spike(s).data,1);
            dev_signal = spike(s).data-median_signal;
            all_dev = cat(3,all_dev,dev_signal);

            xpoints = spike(s).times(2)-spike(s).times(1);%size(spike(s).data,1);
            figure
            set(gcf,'position',[120 183 1440 622])
            plot(linspace(0,spike(s).times(2)-spike(s).times(1),size(spike(s).data,1)),...
                spike(s).data(:,1),'k')
            hold on
            text(xpoints + 0.1,spike(s).data(end,1),spike(s).chLabels{1});
            offset = 0;
            for ich = 2:size(spike(s).data,2)
                offset = offset - (max(spike(s).data(:,ich))-min(spike(s).data(:,ich-1)));
                plot(linspace(0,spike(s).times(2)-spike(s).times(1),size(spike(s).data,1)),...
                    spike(s).data(:,ich)+ offset,'k')
                text(xpoints + 0.1,spike(s).data(end,ich)+offset,spike(s).chLabels{ich});
            end
            xlabel('Time (s)')
            title(sprintf('Time %1.1f s',spike(s).time))
            set(gca,'fontsize',20)
            savefig([plot_folder,ptname,sprintf('_%d',s)]);
            close(gcf)
        end

        avg_dev = mean(all_dev,3);
        figure
        set(gcf,'position',[120 183 1440 622])
        plot(linspace(0,spike(s).times(2)-spike(s).times(1),size(spike(s).data,1)),...
                avg_dev(:,1),'k')
        hold on
        text(xpoints + 0.1,spike(s).data(end,1),spike(s).chLabels{1});
        offset = 0;
        for ich = 2:size(spike(s).data,2)
            offset = offset - (max(avg_dev(:,ich))-min(avg_dev(:,ich-1)));
            plot(linspace(0,spike(s).times(2)-spike(s).times(1),size(spike(s).data,1)),...
                avg_dev(:,ich)+ offset,'k')
            text(xpoints + 0.1,avg_dev(end,ich)+offset,spike(s).chLabels{ich});
        end
        xlabel('Time (s)')
        title(sprintf('Average deviation'))
        set(gca,'fontsize',20)
        savefig([plot_folder,ptname,sprintf('_avg')]);
        close(gcf)
        
    elseif avg == 1
        
        all_z = zeros(size(spike(1).data,1),length(spike));
        
        for s = 1:length(spike)

            mean_signal = mean(spike(s).data,1);
            z_signal = (spike(s).data-mean_signal)./std(spike(s).data);
            z_signal_avg_all_chs = mean(z_signal(:,spike(s).involved),2);
            all_z(:,s) = z_signal_avg_all_chs;
            
        end
        
        all_z_avg = mean(all_z,2);
        figure
        set(gcf,'position',[120 183 1000 300])
        plot(all_z_avg);
        pause
        close(gcf)
        
    end
        
end

end