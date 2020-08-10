function powers = get_power(X,fs,freq_bands)

powers = nan(size(freq_bands,1),1);

% subtract baseline
X = X - median(X);

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

for f = 1:size(freq_bands,1)
    powers(f) = sum(P(freqs>=freq_bands(f,1) & freqs<=freq_bands(f,2)));
end

end