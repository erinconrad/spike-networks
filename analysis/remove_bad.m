function [clean,clean_spike_power,clean_pre] = ...
    remove_bad(bad,orig,spike_power,pre_spike)

clean = orig;
clean_spike_power = spike_power;
clean_pre = pre_spike;

nfreq = length(orig.freq);
np = length(orig.freq(1).pt);

%% Clean main metric
for f = 1:nfreq
    for p = 1:np
        for sp = 1:2
            data = orig.freq(f).pt(p).sp_or_not(sp).data;
            bad_spikes = bad.pt(p).sp_or_not(sp).bad;
            data(bad_spikes,:,:) = [];
            clean.freq(f).pt(p).sp_or_not(sp).data = data;
        end
        
    end
end

%% Clean spike power
for p = 1:np
    bad_spikes = bad.pt(p).sp_or_not(1).bad;
    clean_spike_power(p).spike_powers(bad_spikes,:) = [];
end

%% Clean pre-spike
for p = 1:np
    bad_spikes = bad.pt(p).sp_or_not(1).bad;
    clean_pre(p).before_rise(bad_spikes,:) = [];
end

end