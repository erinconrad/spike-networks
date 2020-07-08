function times = realign_times(nchunks,surround_time)

total_time = surround_time*2;
peak_window = nchunks/2+1;

times = 1:nchunks;
times = times - peak_window; %realign peak to 0
times = times*total_time/nchunks; % adjust for time window

end