function ers_change

%% Parameters

freq_bands = [0 4;... %delta
    4 8;...%theta
    8 12;...% alpha
    12 24;... %beta
    30 40;... % low gamma
    96 106;... % high gamma
    106 256;... %ultra-high
    0 256;... %broadband    
    ]; 
freq_names = {'delta','theta','alpha','beta','low_gamma',...
    'high_gamma','ultra_high','broadband'};

%% Get file locations, load spike times and pt structure
locations = spike_network_files;
main_folder = locations.main_folder;
eeg_folder = [main_folder,'results/eeg_data/'];
results_folder = [main_folder,'results/'];
data_folder = [main_folder,'data/'];
script_folder = locations.script_folder;
addpath(genpath(script_folder));
pt_file = [data_folder,'spike_structures/pt.mat'];
ers_folder = [results_folder,'ers/'];
out_folder = [results_folder,'plots/'];

if exist(out_folder,'dir') == 0
    mkdir(out_folder);
end


% Load spike file for one patient to get the surround time
spike = load([eeg_folder,'HUP074_eeg.mat']);
spike = spike.spike;
surround_time = spike(1).surround_time;


% Loop through time scales
n_freq_abs = 0;
max_F = 0;
network_count = 1;


% Loop through time scales
time_listing = dir(ers_folder);
time_count = 0;

for k = 1:length(time_listing)
    time_name= time_listing(k).name;
    time_window = str2num(time_name);

    % Skip if . or ..
    if strcmp(time_name,'.') == 1 || strcmp(time_name,'..') == 1
        continue
    end

    % Skip if not a directory
    if time_listing(k).isdir == 0, continue; end

    time_count = time_count + 1;
    stats(network_count).time(time_count).name = time_name;
    stats(network_count).time(time_count).time_window = time_window;
    time_folder = [ers_folder,time_name,'/'];

    pt_listing = dir([time_folder,'*.mat']);

    % load one to get nfreq
    ers = load([time_folder,pt_listing(1).name]);
    ers = ers.ers;
    nfreq = size(ers.freq_bands,1);
    if n_freq_abs < nfreq
        n_freq_abs = nfreq;
    end

    for f = 1:nfreq
        stats(network_count).time(time_count).freq(f).avg_ers_involved = ...
            nan(length(pt_listing),surround_time*2/time_window);

        stats(network_count).time(time_count).freq(f).t_all = ...
            nan(length(pt_listing),surround_time*2/time_window);
        
        stats(network_count).time(time_count).freq(f).p_all = ...
            nan(length(pt_listing),surround_time*2/time_window);
    end

    pt_names = {};
    
    % loop through pts
    for i = 1:length(pt_listing)

        pt_name = pt_listing(i).name;
        pt_name_pt = strsplit(pt_name,'_');
        pt_name_pt = pt_name_pt{1};
        pt_names = [pt_names;pt_name_pt];

        % load pt file
        ers = load([time_folder,pt_name]);
        ers = ers.ers;

        for f = 1:nfreq
            stats(network_count).time(time_count).freq(f).name = ers.freq_names{f};
            stats(network_count).time(time_count).freq(f).ers = ers.powers(:,:,f);

            ers_temp = ers.powers(:,:,f);
            
            % Do t-test
            for t = 2:size(ers_temp,2)
                [~,p,~,stats1] = ttest(ers_temp(:,1),ers_temp(:,t));
                stats(network_count).time(time_count).freq(f).t_all(i,t) = stats1.tstat;
                stats(network_count).time(time_count).freq(f).p_all(i,t) = p;
            end
            
            
            
        end


    end
    
    
    % combine t stats across patients to get group-level significance
    for f = 1:nfreq
        stats(network_count).time(time_count).freq(f).group_p = nan(surround_time*2/time_window,1);
        for t = 2:surround_time*2/time_window
            [~,p] = ttest(stats(network_count).time(time_count).freq(f).t_all(:,t));
            stats(network_count).time(time_count).freq(f).group_p(t) = p;
        end
        
        if 0
        fprintf('\nFor %s frequency:\n',freq_names{f});
       (array2table(stats(network_count).time(time_count).freq(f).p_all,...
            'RowNames',pt_names))
        end
    end
    %}
    
    %{
    for f = 1:nfreq
        stats(network_count).time(time_count).freq(f).group_p = nan(surround_time*2/time_window,1);
        for t = 2:surround_time*2/time_window
            
            [~,p] = ttest(stats(network_count).time(time_count).freq(f).avg_ers_involved(:,t),...
                stats(network_count).time(time_count).freq(f).avg_ers_involved(:,1));
         
            stats(network_count).time(time_count).freq(f).group_p(t) = p;
        end
    end
    %}

end
 
%% Plot
% Initialize figure
%{
nfreq columns and 2 (2 time scales) rows
%}
figure
set(gcf,'position',[1 100 1399 500])
[ha, pos] = tight_subplot(time_count, n_freq_abs, [0.04 0.01], [0.13 0.07], [0.06 0.01]);

for t = 1:time_count
    % change times for x axis
    nchunks = size(stats.time(t).freq(1).t_all,2);
    times = realign_times(nchunks,surround_time);
    
    nfreq = length(stats.time(t).freq);
    for f = 1:nfreq
        
        % this adds the number of frequencies + 1 if it's on the 2nd
        % time point (to move down a row), and it adds which frequency
        sp = (n_freq_abs)*(t-1) + f;
        axes(ha(sp));
        
        ers_curr = stats.time(t).freq(f).ers;
        z_curr = (ers_curr-mean(ers_curr,2))./std(ers_curr,0,2);
        
        % loop over patients and plot
        for i = 1:size(ers_curr,1)
            
           
            plot(times,z_curr(i,:),'ko');
            hold on
        end
        
        p = stats.time(t).freq(f).group_p;

        
        % plot mean ers across patients
        for tt = 1:size(ers_curr,2)
            
            text_out = get_asterisks(p(tt),(nchunks-1)*(n_freq_abs));
            
            if strcmp(text_out,'')==1
                plot([times(tt)-0.25 times(tt)+0.25],...
                    [mean(z_curr(:,tt)) mean(z_curr(:,tt))],...
                    'k','linewidth',4);
            else
                plot([times(tt)-0.25 times(tt)+0.25],...
                    [mean(z_curr(:,tt)) mean(z_curr(:,tt))],...
                    'g','linewidth',4);
            end
            
            %{
            if strcmp(text_out,'')==1
                plot([times(tt)-0.25 times(tt)+0.25],...
                    [mean(ers_curr(:,tt)) mean(ers_curr(:,tt))],...
                    'k','linewidth',3);
            else
                plot([times(tt)-0.25 times(tt)+0.25],...
                    [mean(ers_curr(:,tt)) mean(ers_curr(:,tt))],...
                    'g','linewidth',3);
            end
            %}
        end
        
        % asterisks
        %{
        for tt = 2:length(p)
             % should I also adjust by nfreq?
            text(times(tt),max(ers_curr(:,tt))+0.5,sprintf('%s',text_out),'fontsize',20,...
                    'horizontalalignment','center')
        end
        %}
        
        ylim([-2 4])
            
        %{
        if t == 2 && f == 4
             xlabel('Time relative to spike peak (s)')
        end 
        %}
        
        if t == 1
        title(sprintf('%s',...
                strrep(stats.time(t).freq(f).name,'_',' ')))
        end
        
        
        if f == 1 && t == 1
            text(-0.35,-0.15,sprintf('Power in specified\nfrequency band'),...
                'HorizontalAlignment','Center','FontSize',20,'rotation',90,...
                'Units','normalized')
           % ylabel(sprintf('Power in specified\nfrequency band'));
        end
        %}
            
        if t == 1, xticklabels([]); end
        if f~=1, yticklabels([]); end
            
        set(gca,'fontsize',20);

    end

end


annotation('textbox',[0.48 0.03 0.05 0.05],'string',...
    sprintf('Time relative to spike peak (s)'),...
    'HorizontalAlignment','Center','FontSize',20,'linestyle','none',...
    'fitboxtotext','on')



print([out_folder,'ers_change'],gcf,'-depsc');
    

%% Say the patients with significant pre-spike rise
midpoint = nchunks/2;
for t = 1:(time_count)
    fprintf('\n for time window %s:\n\n',time_listing(t).name);
    for i = 1:length(pt_names)
        fprintf('\n%s had significant pre-spike ERS for:',pt_names{i});
        for f = 1:nfreq
            for tt = 1:midpoint - 1
            
                % Get the p-value
                p = stats.time(t).freq(f).p_all(i,tt);
                
                if p < 0.05/nfreq/(nchunks-1)
                    fprintf('\n%s time %d.\n',freq_names{f},tt);
                end

            end
        end
        fprintf('\n\n\n');
    end
end


end