function sac_batch(filenames, pathname, chanlist)
%usage: sac_batch(filenames,pathname,chanlist)
%       filenames is a cell array of strings corresponding to filenames of
%       nev files
%       pathname is a string corresponding to the path the files reside at
%       pathname must have a trailing slash if a directory
%       chanlist is an array of channels to sort (default is 1:96)
%       e.g. sac_batch({'33.nev','35.nev','37.nev'},'/tmp/565r001p')

global FileInfo %contains infrmation about the file/s
global WaveformInfo %Contains Waveform information
global ClusterInfo % Contains channel/classification information
global Handles % Plot/menu handles

FileInfo=struct('filename',[],'format','.nev','HeaderSize',0,...
    'PacketSize',0,'ActiveChannels',[],'PacketOrder',uint8([]),...
    'SpikesNumber',[],'BytesPerSample',0); %initialize FileInfo structure
        
if nargin < 2
    pathname = '';
end

for i=1:length(filenames)
    FileInfo(i).filename=[pathname  filenames{i}] ;
end

ActiveChannelList=[];
for j=1:length(FileInfo)
    sac_nevscan(j);
    ActiveChannelList=union(FileInfo(j).ActiveChannels, ActiveChannelList);
end 

if length(unique([FileInfo(:).PacketSize]))~=1
    error('Variable Packet Sizes');
end

if (nargin < 3) % keep only active channels <= 96 by default
    %chanlist = [1:length(ActiveChannelList)];
    chanlist = ActiveChannelList(find(ActiveChannelList<=96));
end

% only run channels that exist in the file
newchanlist = intersect(chanlist,ActiveChannelList)';

for i=newchanlist
    %for i=1:length(ActiveChannelList)
    %i2=ActiveChannelList(i);
    for j=1:length(FileInfo)
        howmany(i,j)=FileInfo(j).SpikesNumber(i);
    end
    ChannelString{i}=[num2str(i) ' (' num2str(howmany(i,:)) ')'];
end
%set(Handles.channel,'String',ChannelString);
%set(Handles.channel,'UserData',{ActiveChannelList,howmany});
%set(Handles.channel,'Value',1);

for i=newchanlist
%for i=1:length(ActiveChannelList)
    try
        disp(sprintf('\nChannel %i',i));
        %%%%%%%%%%get
        GetNumber=[50 40]; % [Number of Waveform groups, Number of Waveform loaded per group]
	% 50 40 (2000) is default, could try 8000 instead with next line
        %GetNumber=[200 40]; % [Number of Waveform groups, Number of Waveform loaded per group]
        ChannelNumber = i; % ActiveChannelList(i);
        TotalNumber = sum(howmany(i,:));

        disp(sprintf('loading...'));
        sac_get_waveforms_batch(ChannelNumber,TotalNumber,GetNumber);
        unitsmenuString=cat(2,{'noise'},{'unsorted'},num2cell(num2str(setdiff(unique(WaveformInfo.Unit),[0 255]),'%1d')));
    %    set(Handles.unitsmenu,'String',unitsmenuString);
    %    set(Handles.unitsmenu,'Value',2);
        ClusterInfo=[];
        sac_align(0);
        ClusterInfo.Units=unique(WaveformInfo.Unit);
        j=0;
        for j=1:length(ClusterInfo.Units)
            ClusterInfo.Centers(j,:)=mean(WaveformInfo.Waveforms(find(WaveformInfo.Unit==ClusterInfo.Units(j)),:),1);
            ClusterInfo.Sigma{j}=cov(WaveformInfo.Waveforms(find(WaveformInfo.Unit==ClusterInfo.Units(j)),:));
            ClusterInfo.Proportions=length(find(WaveformInfo.Unit==ClusterInfo.Units(j)))/length(WaveformInfo.Unit);
        end
        ClusterInfo.ChannelNumber=ChannelNumber;
        ClusterInfo.microVperBit=FileInfo(1).nVperBit(ChannelNumber)/1000;
        ClusterInfo.Threshold=FileInfo(1).LowThresh(ClusterInfo.ChannelNumber)/FileInfo(1).nVperBit(ClusterInfo.ChannelNumber)*1000;
    %    sac_gui redraw
        disp(sprintf('sorting...'));

        %%%%%%%%%sort        
        unitnum=sac_t_master_batch; %the actual computation
    %    unit_num=max(mod(ClusterInfo.Units,255));
    %    unitsmenuString=cat(2,{'noise'},{'unsorted'},num2cell(num2str(1:unit_num,'%1d')));
    %    set(Handles.unitsmenu,'String',unitsmenuString);
    %    set(Handles.unitsmenu,'Value',min([2,max(unique(WaveformInfo.Unit))+1]));
    %    sac_gui redraw
        disp(sprintf('writing...'));

        %%%%%%%%%write
        sac_write_batch;
    catch
        disp(['Errors on channel ' num2str(i)]);
	disp(['Number of Waveforms = ' num2str(size(WaveformInfo.Waveforms,1))]);
    end
end
