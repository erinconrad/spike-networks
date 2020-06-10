function coherence = erin_alt_coherence(x,y)
    
    gxy = csd(x,y);
    gxx = csd(x,x);
    gyy = csd(y,y);
    
    coherence=  abs(gxy).^2./((gxx).*(gyy));
    
    
    

end


function out = csd(x,y)

if 1
    r = xcorr(x,y);
    
    out= fft(r);
else
    
    
    
    fftx = fft(x);
    ffty = fft(y);
    out = fftx.*conj(ffty);
    
end
    
end