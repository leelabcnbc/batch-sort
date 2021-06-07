function sac_write

global FileInfo;
global ClusterInfo;
global WaveformInfo;

if exist('WaveformInfo') & isfield(WaveformInfo, 'Waveforms')
    store=WaveformInfo.Waveforms;%save Waveforms in temporary variable
else
    store=[];
end

[g,d]=size(ClusterInfo.Centers);
nu=ClusterInfo.nu;
Pj=sparse(ClusterInfo.Proportions);
detSigma=[];
CinvSqrt=[];
for i=1:g
    CinvSqrt{i}=sqrtm(pinv(ClusterInfo.Sigma{i}));
    %detSigma(i)=det(CinvSqrt{i});
    detSigma(i)=1/sqrt(det(ClusterInfo.Sigma{i}));
end
detSigma=diag(sparse(detSigma));
if ~isreal(detSigma)
    warning('Singular covariance')
end
h=waitbar(0,'Writing Files...','Position',[250 300 270 60]);
filenum=length(FileInfo);

%customize according to file type
switch FileInfo(1).format
case '.nev'
    ByteFormat=['int' num2str(FileInfo(1).BytesPerSample*8)];
    UnitFormat='uchar';
    WaveformSize=(FileInfo(1).PacketSize-8)/FileInfo(1).BytesPerSample;% this line changed Feb 5,2003                
    PacketSize=FileInfo(1).PacketSize/FileInfo(1).BytesPerSample;
    Shift1= 8/FileInfo(1).BytesPerSample; %shift to Waveform field (in Packets)
    Shift2= 6;%shift to Unit field (in Bytes)
    UnitLength=1;
    noisecluster=255;
case '.plx'
    ByteFormat='int16';
    UnitFormat='int16';
    WaveformSize=FileInfo(1).PacketSize;
    Shift1=8;%16 bytes to Waveform start
    Shift2=10;
    UnitLength=2;
    ClusterInfo.Units(find(ClusterInfo.Units==255))=0;%plexon programs don't like unit 255
    noisecluster=0;
end

for j=1:filenum %Loop across files
    filename=FileInfo(j).filename;
    HeaderSize=FileInfo(j).HeaderSize;
    fid=fopen(filename,'r+','l');
    Locations=FileInfo(j).Locations;%Locations are in bytes relative to header end
    chan_Locations=Locations(find(FileInfo(j).PacketOrder==ClusterInfo.ChannelNumber));%packets from this channel
    num_to_read=length(chan_Locations);
    Units=zeros(1,num_to_read);
    n=min(3000,round(num_to_read*1e6/Locations(end))); %# of waveforms read at a time
    howmany=ceil(num_to_read/n);%# of read cycles
    current=1;
    
    for i=1:howmany
        
        % Read in the Waveforms
        if i==howmany% last group is smaller
            n=num_to_read-current+1; %line changed August 12, 2003
        end
        waitbar(i/howmany,h);
        relevant_Locations=chan_Locations([current:(current+n-1)])';%Location of packets to read
        fseek(fid,relevant_Locations(1)+HeaderSize,-1);%skip to start location
        relative_Locations=(relevant_Locations-relevant_Locations(1))/FileInfo(j).BytesPerSample;%relative coordinates in sample units
        [chunk,siz]=fread(fid, relative_Locations(end)+Shift1+WaveformSize,['*' ByteFormat]); %read the data in
        Loc_matrix=ones(WaveformSize,1)*relative_Locations+[Shift1+(1:WaveformSize)]'*ones(1,length(relative_Locations));
        WaveformInfo.Waveforms=double(chunk(Loc_matrix))';%create waveform matrix
        
        %   Calculate the classification
        M=zeros(n,g);
        %                 disp(['sorting chunk of ' num2str(n) ' waveforms'])
        sac_align(size(ClusterInfo.Centers,2));
        for k=1:g
            diffs=WaveformInfo.Waveforms-ones(n,1)*ClusterInfo.Centers(k,:);
            M(:,k)=sum(((CinvSqrt{k}*diffs').^2))'; %Mahalanobis distances 
        end
        Prob=exp(-(nu+d)*log(1+M/nu)/2)*detSigma;
        Z=Prob*diag(Pj)./(sum(Prob*diag(Pj),2)*ones(1,g));
        U=(nu+d)./(nu+M);
        [Y,Classification]=max((Z)',[],1);
        Classification=ClusterInfo.Units(Classification);
        Classification(find(max((Z.*U)',[],1)<0.5))=noisecluster; %outliers are assigned to noise cluster
        Units(current:current+n-1)=Classification; %assign final classification
        current=current+n;
    end   
    
    %write the classification back to file
    fseek(fid,HeaderSize+chan_Locations(1)+Shift2,-1);%skip to start location
    skip=[0; diff(chan_Locations)-UnitLength];%first skip is 0 because fwrite first skips then writes
    for i=1:length(skip)
        fwrite(fid,Units(i),UnitFormat,skip(i));
    end
    fclose(fid);
    
end
close(h);
WaveformInfo.Waveforms=store;
