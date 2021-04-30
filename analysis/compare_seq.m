function compare_seq(all,m,f,out_folder)

alpha = 0.05;
np = length(all(1).data.freq(1).pt);

for p = 1:np
    
    %% Get data
    % nsp x ntimes x nchs
    data = all(m).data.freq(f).pt(p).sp_or_not(1).data;
    
    %% Get pre spike info
    pre = all(m).pre;
    
    %% Get sequence info
    
    %% Get max spike power electrode
    sp_power = all(m).data.freq(f).pt(p).spike_powers;
    [~,max_ch] = max(sp_power,[],2);
    
    %% initialize change array
    change = nan(size(data,1),1);
    
    % Loop over spikes
    for s = 1:size(data,1)
        
        %% Get power in the biggest spike power electrode
        main_ch_data(s,:) = data(s,:,max_ch(s));
        current_pre = pre(p).before_rise(s,:);
        
        %% Remove non-pre times
        last_pre = find(current_pre ==0);
        last_pre = last_pre(1) - 1; % the last pre-spike window
        main_ch_data(s,current_pre==0) = nan;
        
        first = main_ch_data(s,1);
        last = main_ch_data(s,last_pre);
        change(s) = last-first;
    end
    
    %% Get sequ
    
end

end