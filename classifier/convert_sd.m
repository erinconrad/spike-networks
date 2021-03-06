function alt_pre_spike = convert_sd(sig_dev,windows,pre_spike,met,do_rand_ch)

    tcount = 0;
    
    % Loop over time scales
    for t = 1:length(sig_dev)
       % sd = sig_dev(t).is_spike(1);
       % sd_not = sig_dev(t).is_spike(2);
        tw = sig_dev(t).is_spike(1).time_window;
        if ~ismember(tw,windows), continue; end
        tcount = tcount + 1;
        
        % Loop over patients
        for p = 1:length(sig_dev(t).is_spike(1).sig_dev)
            
            alt_pre_spike(p).name = pre_spike(p).name;
            if pre_spike(p).windows.which ~= tw, error('what'); end
            alt_pre_spike(p).windows(tcount).which = tw;
            alt_pre_spike(p).windows(tcount).all_windows = pre_spike(p).windows(tcount).all_windows;
            
            alt_pre_spike(p).windows(tcount).manual_before_rise = pre_spike(p).windows(tcount).before_rise;
            
            % Loop over spike vs not spike
            for sp_id = 1:length(sig_dev(t).is_spike)
                
                % if I'm doing not a spike
                if strcmp(sig_dev(t).is_spike(sp_id).name,'not_spike')
                    
                    % picking a random channel
                    if do_rand_ch == 1
                        dev_all = sig_dev(t).is_spike(sp_id).sig_dev(p).dev_windows_rand;
                        
                    else
                        %average across electrodes
                        dev_all = sig_dev(t).is_spike(sp_id).sig_dev(p).dev_windows_auto;
                                   
                    end
                    
                else
                    % doing spike
                    if strcmp(met,'sd_auto')
                        % Get power in automatically measured biggest dev
                        % channel
                        dev_all = sig_dev(t).is_spike(sp_id).sig_dev(p).dev_windows_auto;
                    else
                        dev_all = sig_dev(t).is_spike(sp_id).sig_dev(p).dev_windows;
                    end
                end
                alt_pre_spike(p).windows(tcount).cons_windows = sig_dev(t).is_spike(1).times;
                if sp_id == 1
                    
                    % Fill spike structure
                    alt_pre_spike(p).windows(tcount).dev.spike = dev_all;
                    first = sig_dev(t).is_spike(sp_id).sig_dev(p).dev_windows_first;
                    other = sig_dev(t).is_spike(sp_id).sig_dev(p).dev_windows_other;
                    alt_pre_spike(p).windows(tcount).dev.first = first;
                    alt_pre_spike(p).windows(tcount).dev.other = other;
                    
                else
                    
                    % Fill not a spike structure
                    alt_pre_spike(p).windows(tcount).dev.not = dev_all;
                end
            end
            
        end
    end

end