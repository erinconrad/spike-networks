function clean = remove_bad(orig,bad)

clean = orig;

nfreq = length(orig.freq);
np = length(orig.freq(1).pt);
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


end