function data = download_eeg(dataName,indices,pwname,just_chs)

% This is a tool to return information from a specified iEEG dataset


%% Unchanging parameters
loginname = 'erinconr';

n = 0;


%% Open and get data
% I am putting this whole thing in a try-catch block inside a while loop
% because sometimes the ieeg.org server fails to give data, and this would
% usually crash the program


while n == 0

    try
        
        session = IEEGSession(dataName, loginname, pwname);
        
        if just_chs == 1
            fs = session.data.sampleRate;
            channelLabels = session.data.channelLabels;
        else
            channelLabels = session.data.channelLabels;
            values = session.data.getvalues(indices,':');
        end
    
        
        n = 1;
        
    catch
        fprintf('Failed to retrieve ieeg.org data, trying again...\n'); 
      
        n = 0; 
    end


end


%% Create struct
if just_chs == 1
    data.chLabels = channelLabels;
    data.fs = fs;
else
    data.values = values;
end


session.delete;
clearvars -except data

end