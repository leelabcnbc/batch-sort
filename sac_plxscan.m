function  sac_plxscan(j)
% scans plexon .plx file. Returns into global variable FileInfo variables 
% such as header size and array containing order of locations of packets 
% containing waveforms

global FileInfo;


n = 0;
npw = 0;
ts = 0;
wave = 0;
fid=fopen([ FileInfo(j).filename],'r');
% fid = fopen(filename, 'r');
if(fid == -1)
    disp('cannot open file');
    return
end



% read file header
header = fread(fid, 64, 'int32');
freq = header(35);  % frequency
ndsp = header(36);  % number of dsp channels
nevents = header(37); % number of external events
nslow = header(38);  % number of slow channels
npw = header(39);  % number of points in wave
npr = header(40);  % number of points before threshold
tscounts = fread(fid, [5, 130], 'int32');
wfcounts = fread(fid, [5, 130], 'int32');
evcounts = fread(fid, [1, 512], 'int32');

% read dsp headers
for i=1:ndsp
    name=fread(fid,32,'*char');
    signame=fread(fid,32,'*char');
    channel=fread(fid,1,'int');
    WFrate=fread(fid,1,'int');
    sig=fread(fid,1,'int');
    ref=fread(fid,1,'int');
    FileInfo(j).nVperBit(i)=(6/2^12)/fread(fid,1,'int')*1e6;
    filter=fread(fid,1,'int');
    FileInfo(j).LowThresh(i)=fread(fid,1,'int');
    fseek(fid,928,'cof');
end

% skip other headers
fseek(fid, 296*nevents + 296*nslow, 'cof');
HeaderSize=ftell(fid);

%initialize before scan
record = 1;
keep=int16(zeros(sum(sum(tscounts)),8));
Locations=zeros(sum(sum(tscounts)),1);
loc=0;
chunk_size = 1e6; 
siz=chunk_size;
chunk=int16(zeros(chunk_size,1));
stat=0;

% wave = zeros(npw, 1);
% wf = zeros(npw, 1);

% read data records
while stat==0
    firstloc=loc;
    [chunk,siz]=fread(fid, chunk_size,'*int16');
    while loc<=firstloc+siz-100
        if chunk(loc+1-firstloc)==1
            keep(record,:)=[chunk((loc+1-firstloc):(loc-firstloc+8))]';
            Locations(record)=loc;
            record=record+1;
        end
        loc=loc+8+double(chunk(loc-firstloc+8));        
    end
    stat=feof(fid);
    fseek(fid,loc*2+HeaderSize,'bof');
    
end
fclose(fid);


%Set file information global
FileInfo(1).SamplingRate=freq/1000;%set the sampling rate
FileInfo(1).ThresholdLocation=npr;
FileInfo(j).PacketOrder=uint8(keep(:,5));
FileInfo(j).Locations=Locations*2;%Locations are in units of bytes, relative to end of Header
FileInfo(j).ActiveChannels=double(unique(keep(:,5)));
FileInfo(j).HeaderSize=HeaderSize;
for i1=1:max( FileInfo(j).ActiveChannels) 
    FileInfo(j).SpikesNumber(i1)=sum(FileInfo(j).PacketOrder==i1);   
end
FileInfo(j).PacketSize=npw;
FileInfo(j).BytesPerSample=2; 
