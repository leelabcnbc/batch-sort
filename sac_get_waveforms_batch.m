function sac_get_waveforms_batch(ChannelNumber,TotalNumber,GetNumber)

% This utility loads the waveforms and times from a nev file into WaveformInfo
% The program allows only for one threshold sign, and no A/D saturation
%
% Last modified 1/12/01 Shy Shoham
% Required modification: thresholds should be determined by nev extended header specs
% and not ad-hoc

global FileInfo;
global WaveformInfo;

GetPackets=GetNumber(1);
GetPacketSize=GetNumber(2);
[GetSize,Method]=min([TotalNumber,prod(GetNumber)]);
g=TotalNumber/GetPackets;
ToRead=zeros(TotalNumber,1);
ToRead(find(rem([1:TotalNumber],g)<GetPacketSize))=1; % mark the subset of packets that are to be read

switch FileInfo(1).format
    case '.nev'
        WaveformSize=(FileInfo(1).PacketSize-8)/FileInfo(1).BytesPerSample;
        ByteLength=['int' num2str(FileInfo(1).BytesPerSample*8)];
    case '.plx'
        WaveformSize=FileInfo(1).PacketSize;
        ByteLength='int16';
end
WaveformInfo.Waveforms=zeros(GetSize,WaveformSize);
WaveformInfo.Times=zeros(GetSize,1);
WaveformInfo.Unit=zeros(GetSize,1);
oldnumPackets=0;
oldnumRead=0;
for j=1:length(FileInfo)
    HeaderSize=FileInfo(j).HeaderSize;
    PacketSize=FileInfo(j).PacketSize;
    PacketNumbers=find(FileInfo(j).PacketOrder==ChannelNumber);
    PacketNumber=length(PacketNumbers);
    PacketNumbers=PacketNumbers(find(ToRead((oldnumPackets+1):(oldnumPackets+PacketNumber))));
    fid=fopen(FileInfo(j).filename,'r','l');
    %     jump=[PacketSize*(diff(PacketNumbers)-1); 0];% array specifying size of jumps
%     jump=diff(FileInfo(j).Locations(PacketNumbers));
    ReadNumber=length(PacketNumbers);
    if ReadNumber>0
        %         fseek(fid,HeaderSize+(PacketNumbers(1)-1)*PacketSize,-1);
        fseek(fid,HeaderSize+FileInfo(j).Locations(PacketNumbers(1)),-1);%move to first relevant packet
        for i=oldnumRead+[1:ReadNumber] 
            switch FileInfo(1).format
                case '.nev'
                    WaveformInfo.Times(i,1)=fread(fid,1,'uint32');
                    PacketIdentifier=fread(fid,1,'uint16');
                    WaveformInfo.Unit(i,1)=fread(fid,1,'uchar');
                    Temp=fread(fid,1,'uchar');
                    WaveformInfo.Waveforms(i,:)=fread(fid,WaveformSize,ByteLength)';
                case '.plx'
                    fseek(fid,4,0);
                    WaveformInfo.Times(i,1)=fread(fid,1,'int32');
                    PacketIdentifier=fread(fid,1,'int16');
                    WaveformInfo.Unit(i,1)=fread(fid,1,'int16');
                    fseek(fid,4,0);
                    WaveformInfo.Waveforms(i,:)=fread(fid,FileInfo(j).PacketSize,ByteLength)';
            end
            fseek(fid,HeaderSize+FileInfo(j).Locations(PacketNumbers(i-oldnumRead)),-1);
%             fseek(fid,jump(i-oldnumRead),'cof');
        end
    end
    oldnumRead=oldnumRead+ReadNumber;
    oldnumPackets=oldnumPackets+PacketNumber;
    fclose(fid);
end   

%Rejecting waveforms that reached A/D saturating
if FileInfo(1).BytesPerSample==1,
    maxval=127;
    minval=-127;
else
    maxval=2047;
    minval=-2048;
end
%threshold=mean(WaveformInfo.Waveforms(:,11));
%thresholdsign=sign(threshold);
%goodIndices=find(sign(WaveformInfo.Waveforms(:,11))==thresholdsign...
%   & max(WaveformInfo.Waveforms')'<maxval & min(WaveformInfo.Waveforms')'> minval );
goodIndices=find(max(WaveformInfo.Waveforms')'<maxval & min(WaveformInfo.Waveforms')'> minval );
WaveformInfo.Times=WaveformInfo.Times(goodIndices);
WaveformInfo.Unit=WaveformInfo.Unit(goodIndices);
WaveformInfo.Unit(WaveformInfo.Unit>8 & WaveformInfo.Unit<255)=0;
WaveformInfo.Waveforms=WaveformInfo.Waveforms(goodIndices,:);