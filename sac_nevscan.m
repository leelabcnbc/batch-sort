function nevscan(j)
global FileInfo;
fid=fopen([ FileInfo(j).filename],'r','l');

%read general header
FileType=fread(fid,8,'uchar');
Version=fread(fid,2,'uchar');
FileFormatAdditional=fread(fid,2,'uchar');
HeaderSize=fread(fid,1,'uint32');
PacketSize=fread(fid,1,'uint32');
TimeReslutionTimeStamps=fread(fid,1,'uint32');
TimeReslutionSamples=fread(fid,1,'uint32');
%unsure about actualtype of TimeOrigin
TimeOrigin=fread(fid,8,'uint16');
Application=fread(fid,32,'uchar');
Comment=fread(fid,256,'uchar');
ExtendedHeaderNumber=fread(fid,1,'ulong');

%read extended headers
for i=1:ExtendedHeaderNumber
    Identifier=char(fread(fid,8,'char'))';
    %modify this later
    switch Identifier
        case 'NEUEVWAV'
            ElecID=fread(fid,1,'uint16');
            PhysConnect=fread(fid,1,'uchar');
            PhysConnectPin=fread(fid,1,'uchar');
            FileInfo(j).nVperBit(ElecID)=fread(fid,1,'uint16');
            EnergyThresh=fread(fid,1,'uint16');
            FileInfo(j).HighThresh(ElecID)=fread(fid,1,'int16');
            FileInfo(j).LowThresh(ElecID)=fread(fid,1,'int16');
            SortedUnits=fread(fid,1,'uchar');
            BytesPerSample=((fread(fid,1,'uchar'))>1)+1;
            temp=fread(fid,10,'uchar');
        otherwise, % added26/7/05 after identifying bug in reading etended headers
            temp=fread(fid,24,'uchar');
    end
end

% Calculate number of packets
fseek(fid,0,1);
FileSize=ftell(fid);
PacketNumber=(FileSize-HeaderSize)/PacketSize;

%initialization
fseek(fid,HeaderSize,-1);
fread(fid,1,'uint32');


%read the data
FileInfo(j).PacketOrder=uint8(fread(fid,PacketNumber,'uint16',PacketSize-2));%read the channel identifiers
fseek(fid,HeaderSize,-1);
Times=fread(fid,PacketNumber,'uint32',PacketSize-4);%read the packet timestamps
fclose(fid);
FileInfo(j).Locations=[0:PacketSize:PacketSize*(PacketNumber-1)]';%The location of packets.
% This line was changed, August 12, 2003

%Determine active channels and number of spikes on each

FileInfo(j).SpikesNumber=zeros(1,255);

for i1=1:double(max(FileInfo(j).PacketOrder))
    FileInfo(j).SpikesNumber(i1)=length(find(FileInfo(j).PacketOrder==i1));
end
FileInfo(j).ActiveChannels=find(FileInfo(j).SpikesNumber);

% The next section deals with finding and purging amplifier pops

ShockThreshold=round(0.9*length(FileInfo(j).ActiveChannels));

if ShockThreshold>30
    suspects=find((Times(ShockThreshold:end)-Times(1:end-ShockThreshold+1))<=2);
    for i2=suspects'
        ShockPackets=find(abs(Times-Times(i2))<2);
        FileInfo(j).PacketOrder(ShockPackets)=0;
    end
    for i1=1:max( FileInfo(j).ActiveChannels)
        FileInfo(j).SpikesNumber(i1)=sum(FileInfo(j).PacketOrder==i1);
    end
end

%Set file information global
FileInfo(1).SamplingRate=30;
FileInfo(1).ThresholdLocation=11;
FileInfo(j).ActiveChannels=find(FileInfo(j).SpikesNumber);
FileInfo(j).HeaderSize=HeaderSize;
FileInfo(j).PacketSize=PacketSize;
FileInfo(j).BytesPerSample=BytesPerSample; %This can have electrode dependent values, but here only set once