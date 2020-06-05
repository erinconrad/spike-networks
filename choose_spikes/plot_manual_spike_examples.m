function plot_manual_spike_examples(whichInListing,avg,just_involved,just_save)

%{
This function plots examples of the EEG data surrounding manually detected
spikes. It will either plot all of the channels or the average signal
z-score of the involved channels.
%}

%% Parameters
num_to_plot = 50; % number of spikes to plot per pt

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


if isempty(whichInListing) == 1
    whichInListing = 1:length(listing);
end

for i = whichInListing
    
    split_name = strsplit(listing(i).name,'_');
    ptname = split_name{1};
    fprintf('Doing %s\n',ptname);
    
    eeg = load([eeg_folder,listing(i).name]);
    spike = eeg.spike;
    
    % Pick a random num_to_plot spikes to plot
    s_to_plot = randsample(1:length(spike),num_to_plot);
    
    all_dev = [];
    fs = spike(1).fs;
    
    if avg == 0
    
        for s = s_to_plot
            which_chs = 1:size(spike(s).data,2);
            if just_involved == 1
                which_chs = which_chs(spike(s).involved);
            end
            
            signal = spike(s).data;
            signal = pre_processing(signal,1,0,1,fs);
            mean_signal = mean(signal,1);
            z_signal = (signal-mean_signal)./std(signal);
            all_dev = cat(3,all_dev,z_signal);

            xpoints = spike(s).times(2)-spike(s).times(1);%size(spike(s).data,1);
            figure
            set(gcf,'position',[120 183 1440 622])
            offset = 0;
            for ich = 1:length(which_chs)
                if ich > 1
                    offset = offset - (max(signal(:,which_chs(ich)))-min(signal(:,which_chs(ich-1))));
                else
                    offset = 0;
                end
                plot(linspace(0,spike(s).times(2)-spike(s).times(1),size(signal,1)),...
                    signal(:,which_chs(ich))+ offset,'k')
                hold on
                text(xpoints + 0.1,signal(end,which_chs(ich))+offset,spike(s).chLabels{which_chs(ich)});
            end
            xlabel('Time (s)')
            title(sprintf('Time %1.1f s for %s',spike(s).time,ptname))
            set(gca,'fontsize',20)
            pause
            %savefig([plot_folder,ptname,sprintf('_%d',s)]);
            close(gcf)
        end
        %{
        avg_dev = mean(all_dev,3);
        figure
        set(gcf,'position',[120 183 1440 622])
        offset = 0;
        for ich = 1:length(which_chs)
            if ich > 1
                offset = offset - (max(avg_dev(:,which_chs(ich)))-min(avg_dev(:,which_chs(ich-1))));
            else
                offset = 0;
            end
            plot(linspace(0,spike(s).times(2)-spike(s).times(1),size(spike(s).data,1)),...
                avg_dev(:,ich)+ offset,'k')
            text(xpoints + 0.1,avg_dev(end,which_chs(ich))+offset,spike(s).chLabels{which_chs(ich)});
        end
        xlabel('Time (s)')
        title(sprintf('Average deviation'))
        set(gca,'fontsize',20)
        pause
        %savefig([plot_folder,ptname,sprintf('_avg')]);
        close(gcf)
        %}
        
    elseif avg == 1
        
        all_z = zeros(size(spike(1).data,1),length(spike));
        
        for s = 1:length(spike)
            signal = spike(s).data;
            signal = pre_processing(signal,1,0,1,fs);
            mean_signal = mean(signal,1);
            z_signal = abs(signal-mean_signal)./std(signal);
            z_signal_avg_all_chs = mean(z_signal(:,spike(s).involved),2);
            all_z(:,s) = z_signal_avg_all_chs;
            
        end
        
        all_z_avg = nanmean(all_z,2);
        figure
        set(gcf,'position',[120 183 1000 300])
        plot([1:size(all_z_avg,1)]/fs,all_z_avg);
        xlabel('Time (s)')
        ylabel('Z score')
        title(sprintf('Average signal deviation surrounding spikes for %s\n(involved channels only)',ptname))
        xlim([0 6])
        set(gca,'fontsize',20)
        if just_save == 0
            pause
        end
        print(gcf,[plot_folder,ptname,sprintf('avg')],'-depsc')
        close(gcf)
        
    end
        
end

end