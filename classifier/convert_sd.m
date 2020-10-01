function alt_pre_spike = convert_sd(sig_dev,windows,pre_spike)

    tcount = 0;
    for t = 1:length(sig_dev)
        sd = sig_dev(t).is_spike(1);
        sd_not = sig_dev(t).is_spike(2);
        tw = sd.time_window;
        if ~ismember(tw,windows), continue; end
        tcount = tcount + 1;
        for p = 1:length(sd.sig_dev)
            
            
            pval = sd.sig_dev(p).p;
            h = pval < 0.1;
            %}
            
            sp_dev = sd.sig_dev(p).dev_windows;
            thresh = 1;
            mean_sp_dev = mean(sp_dev,1);
            std_sp_dev = std(sp_dev,1);
            
            
            
            alt_pre_spike(p).name = pre_spike(p).name;
            if pre_spike(p).windows.which ~= tw, error('what'); end
            alt_pre_spike(p).windows(tcount).which = tw;
            alt_pre_spike(p).windows(tcount).all_windows = pre_spike(p).windows(tcount).all_windows;
            alt_pre_spike(p).windows(tcount).before_rise = ones(size(pre_spike(p).windows(tcount).before_rise));
            
            alt_pre_spike(p).windows(tcount).before_rise(:,h==1) = ...
                0;
            alt_pre_spike(p).windows(tcount).before_rise(:,1) = 1;
            alt_pre_spike(p).windows(tcount).dev.spike = sd.sig_dev(p).dev_windows;
            alt_pre_spike(p).windows(tcount).dev.not = sd_not.sig_dev(p).dev_windows;
            alt_pre_spike(p).windows(tcount).manual_before_rise = pre_spike(p).windows(tcount).before_rise;
            
            alt_pre_spike(p).windows(tcount).cons = ones(size(pre_spike(p).windows(tcount).before_rise));
            alt_pre_spike(p).windows(tcount).cons(:,(sig_dev(t).tests.unpaired.sig==1)') = 0;
            alt_pre_spike(p).windows(tcount).cons(:,1) = 1;
        end
    end

end