function powers = get_power(X,fs,freq_bands)

powers = nan(size(freq_bands,1),1);

% subtract baseline (I do this earlier at the level of the full spike time
% period rather than individual windows)
%X = X - median(X);

%% fft approach

% Calculate fft
Y = fft(X);

% Get power
P = abs(Y).^2;
freqs = linspace(0,fs,length(P)+1);
freqs = freqs(1:end-1);

% Take first half
P = P(1:ceil(length(P)/2));
freqs = freqs(1:ceil(length(freqs)/2));

% Perform linear interpolation to smooth this out
freqs = [freqs,256];
P = [P;P(end)];
finterp = [0:2:256];
Pinterp = interp1(freqs',P,finterp);


for f = 1:size(freq_bands,1)
    powers(f) = sum(Pinterp(finterp>=freq_bands(f,1) & finterp<=freq_bands(f,2)));
end

end