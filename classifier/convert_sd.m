function alt_pre_spike = convert_sd(sig_dev,windows,pre_spike)

    
    for t = 1:length(sig_dev)
        sd = sig_dev(t).is_spike(1);
        tw = sd.time_window;
        if ~ismember(tw,windows), continue; end
        for p = 1:length(sd.sig_dev)
            
            pval = sd.sig_dev(p).p;
            h = pval < 0.05;
            
            alt_pre_spike(p).name = pre_spike(p).name;
            if pre_spike(p).windows.which ~= tw, error('what'); end
            alt_pre_spike(p).windows(t).which = tw;
            alt_pre_spike(p).windows(t).all_windows = pre_spike(p).windows(t).all_windows;
            alt_pre_spike(p).windows(t).before_rise = ones(size(pre_spike(p).windows(t).before_rise));
            
            alt_pre_spike(p).windows(t).before_rise(:,h==1) = ...
                0;
            alt_pre_spike(p).windows(t).before_rise(:,1) = 1;
        end
    end

end