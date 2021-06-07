function write_nev(Waveforms,Times,Freq,nV_per_bit)

% write_nev(Waveforms,Times,Freq,nV_per_bit)
%
% This function writes a nev file containing a single channel 
% Waveforms - Waveforms to be written (Matrix: # of Waveofrms x Waveform length)
% Times - vector containing timestamps (one for each waveform)
% Freq - sampling frequency (Hz) used (Timestamps should be in units of 1/Freq).
% nV_per_bit - gain used. Waveforms*nV_per_bit should gives the actual amplitude in nanoVolts
% 
% If you have a data file open in SAC and you want to write one channel
% into a file, use the following two commands:
% > global ClusterInfo WaveformInfo FileInfo
% > write_nev(WaveformInfo.Waveforms,WaveformInfo.Times,FileInfo.SamplingRate,FileInfo.nVperBit(ClusterInfo.ChannelNumber));


[count,waveform_length]=size(Waveforms);

[filename,pathname]=uiputfile('*.nev','name of file to write');
fid=fopen([pathname filename],'w','l');
fwrite(fid,'NEURALEV','char');%FileType,8
fwrite(fid,[2 0],'uchar');%Version,2
fwrite(fid,[0 0],'char');%format additional,2
fwrite(fid,368,'uint32');% Header Size,uint32. 336+32*#of extended headers
fwrite(fid,waveform_length*2+8,'uint32');% Packet Size, uint32
fwrite(fid,Freq,'uint32');%sampling frequency of time stamps, uint32
fwrite(fid,Freq,'uint32');%sampling frequency of samples
fwrite(fid,zeros(1,8),'uint16');%Time origin,8xuint16
fwrite(fid,zeros(1,32),'char');%application,32
fwrite(fid,zeros(1,256),'char');%comment,256
fwrite(fid,1,'uint32');%extended header number, uint32

%Write extended header
fwrite(fid,'NEUEVWAV','char'); %Packet type identifier, 8
fwrite(fid,1,'uint16'); % Electrode ID  
fwrite(fid,1,'uchar'); % Physical connector
fwrite(fid,1,'uchar'); % Physical connector pin
fwrite(fid,nV_per_bit,'uint16'); % Gain, nV per bit
fwrite(fid,0,'uint16'); % EnergyThresh
fwrite(fid,0,'int16'); % HighThresh
fwrite(fid,0,'int16'); % LowThresh
fwrite(fid,0,'uchar'); % SortedUnits
fwrite(fid,2,'uchar'); % BytesPerSample
fwrite(fid,zeros(10,1),'uchar');% ?

%Write data packets
h=waitbar(0,'writing file...');
for i1=1:count
   if mod(i1,100)==0
      waitbar(i1/count)
   end
   fwrite(fid,Times(i1),'uint32');%Timestamp  
   fwrite(fid,1,'uint16');%Channel #
   fwrite(fid,0,'uchar');%Unit #
   fwrite(fid,0,'uchar');%Temporary place holder
   fwrite(fid,Waveforms(i1,:),'int16');%Waveform
end
close(h);
fclose(fid);

