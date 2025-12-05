function MUAartrem(path)
% Note that this code requires the installation of the TDT MATLAB offline
% analysis tools: https://www.tdt.com/docs/sdk/offline-data-analysis/offline-data-matlab/
% as well as the ARtACS software, DOI: 10.5281/zenodo.596319
% https://github.com/agricolab/ARtACS

format long

%% Set up folder/file selection
selpath = path;

data = TDTbin2mat(selpath);


%% Extract needed data

Fs_MUA = double(data.streams.Wav2.fs); %sampling frequency for MUA channel
Fs_stim = double(data.streams.uS1r.fs); %sampling frequency of stimulation channel

MUAdata = double(data.streams.Wav2.data); %MUA data
stimdata = [data.streams.MonA.data(1,:); data.streams.MonA.data(2,:)]; %have to check on TDT system whether it is meant to be addition or subtraction and what is the unit..?

a = size(data.epocs);
if a(1) == 1
    disp('yes')

    %1-min sACS paradigms have a size of 3, 20-min paradigms have a size
    %of 1
    b = size(data.epocs.Wfrq.onset);
    if b(1) > 1
        freq = double(int32(data.epocs.Wfrq.data(2))); %frequency of stimulation
        amp = double(int32(data.epocs.Wap_.data(2)));  %amplitude of stimulation

        stimon = int32(data.epocs.Wfrq.onset(2)*Fs_MUA); %index where stimulation turns on in LFP channel
        stimoff = int32(data.epocs.Wfrq.offset(2)*Fs_MUA); %index where stimulation turns off in LFP channel
    else
        freq = double(int32(data.epocs.Wfrq.data(1))); %frequency of stimulation
        amp = double(int32(data.epocs.Wap_.data(1)));  %amplitude of stimulation

        stimon = int32(data.epocs.Wfrq.onset(1)*Fs_MUA); %index where stimulation turns on in LFP channel
        stimoff = int32(data.epocs.Wfrq.offset(1)*Fs_MUA); %index where


        %% Filtering MUA recording with comb filter and plotting before and after
    MUAdur = MUAdata(:,stimon:stimoff);

    numperiods = 10;


    % Comb filter to remove artifact harmonics
    if freq > 100
        adjustfilt = 10;
    elseif (freq < 100) && (freq > 10)
        adjustfilt = 4;
    elseif freq <= 10
        adjustfilt = 1;
    end

    for i = 1:size(MUAdur,1)
        MUAfilt(i,:) = artacs.dft.causal(MUAdur(i,:),[freq-adjustfilt:freq+adjustfilt],round(Fs_MUA),numperiods);
    end

    

    for i = 1:size(MUAdata,1)
        MUAdata(i,stimon:stimoff) = MUAfilt(i,:);
    end

    %% Saving data for next step in processing (Spyking -> Python)
    MUAsave.MUAdata = MUAdata;
    MUAsave.Fs = Fs_MUA;
    MUAsave.amp = amp;
    MUAsave.freq = freq;
    MUAsave.stimon = [stimon stimoff];
    MUAsave.stimdata = stimdata;
    save(append(selpath, '/MUAfiltdata.mat'), 'MUAsave', '-v7.3')

else
    disp('no')
    MUAsave.MUAdata = MUAdata;
    MUAsave.Fs = Fs_MUA;
    save(append(selpath, '/MUAfiltdata.mat'), 'MUAsave', '-v7.3')
end

end
