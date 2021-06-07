function sac_gui(action)
global FileInfo %contains infrmation about the file/s
global WaveformInfo %Contains Waveform information
global ClusterInfo % Contains channel/classification information
global Handles % Plot/menu handles


colors=[0 0 0;0.3 0.3 0.3;1 0 0; 0 1 0; 0 0 1; 0.5 0.5 0; 0.5 0 0.5;0 0.7 0.7;1 1 0; 1 0 1;]; 
colornames={'Black','Silver','Red','Green','Blue','Bronze','Magenta','Cyan','Yellow', 'Pink'};
Handles.channel=findobj(gcbf,'Tag','ChannelNo');
Handles.waveforms=findobj(gcbf,'Tag','waveforms');
Handles.times=findobj(gcbf,'Tag','times');
Handles.PCA=findobj(gcbf,'Tag','PCA');
Handles.PCA_num=findobj(gcbf,'Tag','PC_num');
Handles.Penalty=findobj(gcbf,'Tag','Penalty');
Handles.unitsmenu=findobj(gcbf,'Tag','unitmenu');
if ~isempty(WaveformInfo)
    unitmap=[{255} {0} num2cell(setdiff(unique(WaveformInfo.Unit),[0 255]))];
end

switch(action)
    
    case 'openfile'
        
        [filename,pathname]=uigetfile({'*.nev;*.plx'}, 'Choose file');
        if filename==0
            return;
        end
        [temp, extension]=strtok(filename,'.'); %determine data type
        FileInfo=struct('filename',[pathname filename],'format',extension, 'HeaderSize',0,...
            'PacketSize',0,'ActiveChannels',[],'PacketOrder',uint8([]),...
            'SpikesNumber',[]); %initialize FileInfo structure
        
        sac_gui newfile;
        
    case 'multiple' %open multiple files
        
        presentdir=pwd;
        [filenames,pathname]=uigetfiles('*.*', 'Choose files');
        filenumber=length(filenames);
        
        if pathname==0
            return;
        end
        
        FileInfo=FileInfo([]);
        [temp, extension]=strtok(filenames{1},'.'); %determine data type
        FileInfo=struct('filename',[],'format',extension,'HeaderSize',0,...
            'PacketSize',0,'ActiveChannels',[],'PacketOrder',uint8([]),...
            'SpikesNumber',[],'BytesPerSample',0); %initialize FileInfo structure
        
        for i=1:filenumber
            FileInfo(i).filename=[pathname  filenames{i}] ;
        end
        sac_gui newfile
        
    case 'newfile' % open a new file
        
        h=waitbar(0,'Scanning Files...');
        ActiveChannelList=[];
        for j=1:length(FileInfo)
            waitbar((j-1)/length(FileInfo));
            switch FileInfo(1).format
                case '.nev'
                    sac_nevscan(j);
                case '.plx'
                    sac_plxscan(j);
            end
            ActiveChannelList=union(FileInfo(j).ActiveChannels, ActiveChannelList);         
        end 
        waitbar(1);
        close(h);
        
        if length(unique([FileInfo(:).PacketSize]))~=1
            error('Variable Packet Sizes');
        end
        
        for i=1:length(ActiveChannelList)
            i2=ActiveChannelList(i);
            for j=1:length(FileInfo)
                howmany(i,j)=FileInfo(j).SpikesNumber(i2);
            end      
            ChannelString{i}=[num2str(i2) ' (' num2str(howmany(i,:)) ')'];            
        end   
        set(Handles.channel,'String',ChannelString);
        set(Handles.channel,'UserData',{ActiveChannelList,howmany});
        set(Handles.channel,'Value',1);
        
    case 'write'  %Classify all waveforms in this channel
        sac_write;
        
    case 'save_params' % write cluster information into a parameter file
        
        [filename,pathname]=uiputfile('*_ruleset.mat','file to save parameters in');
        present=pwd;
        cd (pathname)
        if exist(filename,'file')~=0
            load([pathname filename],'ruleset');
        end
        cd(present);
        ruleset(ClusterInfo.ChannelNumber).Centers=ClusterInfo.Centers;  
        ruleset(ClusterInfo.ChannelNumber).Sigma=ClusterInfo.Sigma;
        ruleset(ClusterInfo.ChannelNumber).Units=ClusterInfo.Units;
        ruleset(ClusterInfo.ChannelNumber).Proportions=ClusterInfo.Proportions;
        ruleset(ClusterInfo.ChannelNumber).nu=ClusterInfo.nu;
        ruleset(ClusterInfo.ChannelNumber).Version=2;
        save([pathname filename],'ruleset');
        
    case 'apply_params'
        sac_apply_params;
        
        
    case 'get', %getting waveforms from file
        
        GetNumber=[50 40]; % [Number of Waveform groups, Number of Waveform loaded per group]
        ChannelInfo= get(Handles.channel,'UserData');
        ChannelNumber=ChannelInfo{1}(get(Handles.channel,'Value'));
        TotalNumber=sum(ChannelInfo{2}(get(Handles.channel,'Value'),:));
        sac_get_waveforms(ChannelNumber,TotalNumber,GetNumber);
        unitsmenuString=cat(2,{'noise'},{'unsorted'},num2cell(num2str(setdiff(unique(WaveformInfo.Unit),[0 255]),'%1d')));
        set(Handles.unitsmenu,'String',unitsmenuString);
        set(Handles.unitsmenu,'Value',2);
        ClusterInfo=[];
        sac_align(0);
        ClusterInfo.Units=unique(WaveformInfo.Unit);
        i=0;
        for i=1:length(ClusterInfo.Units)
            ClusterInfo.Centers(i,:)=mean(WaveformInfo.Waveforms(find(WaveformInfo.Unit==ClusterInfo.Units(i)),:),1);
            ClusterInfo.Sigma{i}=cov(WaveformInfo.Waveforms(find(WaveformInfo.Unit==ClusterInfo.Units(i)),:));
            ClusterInfo.Proportions=length(find(WaveformInfo.Unit==ClusterInfo.Units(i)))/length(WaveformInfo.Unit);
        end
        ClusterInfo.ChannelNumber=ChannelNumber;
        ClusterInfo.microVperBit=FileInfo(1).nVperBit(ChannelNumber)/1000;
        ClusterInfo.Threshold=FileInfo(1).LowThresh(ClusterInfo.ChannelNumber)/FileInfo(1).nVperBit(ClusterInfo.ChannelNumber)*1000;
        sac_gui redraw
        
        
        
    case 'sort',
        
        unitnum=sac_t_master; %the actual computation
        unit_num=max(mod(ClusterInfo.Units,255));
        unitsmenuString=cat(2,{'noise'},{'unsorted'},num2cell(num2str(1:unit_num,'%1d')));
        set(Handles.unitsmenu,'String',unitsmenuString);
        set(Handles.unitsmenu,'Value',min([2,max(unique(WaveformInfo.Unit))+1]));
        sac_gui redraw
        
    case 'changemap'
        unit=get(Handles.unitsmenu,'Value');
        unitnumber=length(unitmap);
        colornames_temp=colornames([-1 0 unitmap{3:unitnumber}]+2);
        colornames_temp{unit}='----';
        [selection,ok]=listdlg('PromptString',['Add to ' colornames{unit}],...
            'InitialValue',min(unitnumber,unit+1),'ListString',colornames_temp);
        if ok==1 %add selected units to a new unit and update the unit map
            unitlist=cat(2,unitmap{selection});%list of selected units
            % update the classifications
            WaveformInfo.Unit(find(ismember(WaveformInfo.Unit,unitlist)))=unitmap{unit};    
            % update the cluster->unit map
            ClusterInfo.Units(find(ismember(ClusterInfo.Units,unitlist)))=unitmap{unit};
            %update the unit selection box
            unitsmenuString=cat(2,{'noise'},{'unsorted'},num2cell(num2str(setdiff(unique(WaveformInfo.Unit),[0 255]),'%1d')));
            set(Handles.unitsmenu,'String',unitsmenuString);
            set(Handles.unitsmenu,'Value',2);
            sac_gui redraw    
        end
        
    case 'limits'
        while 1
            m=mean(WaveformInfo.Waveforms);
            [x,y,button]=ginput(1);
            x=x*FileInfo(1).SamplingRate;
            y=y/ClusterInfo.microVperBit;
            if button>1
                break
            end
            smaller=WaveformInfo.Waveforms(:,round(x))<y;
            good=find(abs(smaller-(y<m(round(x)))));
            WaveformInfo.Waveforms=WaveformInfo.Waveforms(good,:);
            WaveformInfo.Unit=WaveformInfo.Unit(good);
            WaveformInfo.Times=WaveformInfo.Times(good);
            sac_gui redraw
        end
        
        
    case 'de-noise'
        ToKeep=find(WaveformInfo.Unit~=255);
        WaveformInfo.Waveforms=WaveformInfo.Waveforms(ToKeep,:);
        WaveformInfo.Times=WaveformInfo.Times(ToKeep);
        WaveformInfo.Unit=WaveformInfo.Unit(ToKeep);
        sac_gui redraw
        
    case 'add-new'
        i=get(Handles.unitsmenu,'Value');
        unit=unitmap{i};
        new_unit=unitmap{end}+1;
        WaveformInfo.Unit(find(WaveformInfo.Unit==unit))=new_unit;
        ClusterInfo.Units(find(ClusterInfo.Units==unit))=new_unit;
        unitsmenuString=cat(2,{'noise'},{'unsorted'},num2cell(num2str(setdiff(unique(WaveformInfo.Unit),[0 255]),'%1d')));
        set(Handles.unitsmenu,'String',unitsmenuString);
        set(Handles.unitsmenu,'Value',length(unitsmenuString));
        sac_gui redraw
        
    case 'redraw'
        unit=get(Handles.unitsmenu,'Value');
        unitnumber=length(unitmap);
        set(gcbf,'CurrentAxes',Handles.waveforms);
        cla;
        axis tight
        set(gca,'Nextplot','add','FontSize',13);
        xlabel('Time(msec)')
        ylabel('microVolts')
        %Waveform plot
        for i=[setdiff(1:unitnumber,unit) unit] %order assures that current unit is drawn last
            relevant=find(ismember(WaveformInfo.Unit,unitmap{i}));    
            
            relevant=find(ismember(WaveformInfo.Unit,unitmap{i}));
            %             if(~isempty(relevant))
            %                 temporary_variable = round(linspace(1, length(relevant) - 2, 1e3));
            %                 relevant = relevant(temporary_variable);
            %             end
            
            if length(relevant)>0
                wh=plot([1:size(WaveformInfo.Waveforms,2)]'*ones(1,length(relevant))/FileInfo(1).SamplingRate, ClusterInfo.microVperBit*WaveformInfo.Waveforms(relevant,:)','Color',colors(mod(unitmap{i}+2,256),:));
            end
        end
        axis tight
        %ISI plot
        set(gcbf,'CurrentAxes',Handles.times);
        cla;
        title(' ')
        set(gca,'Nextplot','replacechildren','FontSize',13);   
        hold on
        shading faceted
        %for i=[setdiff(1:unitnumber,unit) unit] %order assures that current unit is drawn last
        for i=unit
            relevant=find(ismember(WaveformInfo.Unit,unitmap{i}));    
            diffs=diff(WaveformInfo.Times(relevant))/FileInfo(1).SamplingRate;
            if length(diffs((diffs<1e2)&(diffs>0)))>0
                [N,bins]=histc(diffs((diffs<1e2)&(diffs>0)),[0:1:100]);  
                edges=[0:1:100];
                h=bar(edges,N/(length(relevant)),'histc'); 
                axis tight;
                set(h,'EdgeColor',colors(mod(unitmap{i}+2,256),:),'FaceColor','none')
                xlabel('msec')
                ylabel('Pr(ISI)')
            end      
        end
        title([num2str(length(relevant)) ' Waveforms']);
        
        %PCA plot
        set(gcbf,'CurrentAxes',Handles.PCA);
        cla
        axis auto
        set(gca,'Units','characters','FontUnits','points','FontSize',13);
        [pc,score,latent] = princomp(WaveformInfo.Waveforms);
        for i=[setdiff(1:unitnumber,unit) unit] %order assures that current unit is drawn last
            relevant=find(ismember(WaveformInfo.Unit,unitmap{i}));    
            if length(relevant)>0
                wh=plot(score(relevant,1),score(relevant,2), '.','MarkerSize',4, 'Color',colors(mod(unitmap{i}+2,256),:));
            end
        end
        zero=(zeros(size(mean(WaveformInfo.Waveforms)))-mean(WaveformInfo.Waveforms))*pc(:,1:2);
        plot(zero(1),zero(2),'k*','MarkerSize',8)
        xlabel('PC 1');
        ylabel('PC 2');
        if  isfield(ClusterInfo,'Centers')
            clusternum=size(ClusterInfo.Centers,1);
            hold on     
            for i=1:clusternum
                sig1=pc(:,1:2)'*(ClusterInfo.Sigma{i}*pc(:,1:2));%rotate covariance to pc dimensions
                [V,D]=eig(sig1);
                mu=(ClusterInfo.Centers(i,:)-mean(WaveformInfo.Waveforms))*pc(:,1:2);
                ph1(i)=sac_ellipse(2*sqrt(D(1,1)),2*sqrt(D(2,2)),angle(V(1,1)+1i*V(2,1)),mu(1),mu(2));
                set(ph1(i),'LineWidth',2,'Color',colors(mod(ClusterInfo.Units(i)+2,256),:));
            end
        end
        
    case 'change_plot'
        
        set(gcbf,'CurrentAxes',Handles.waveforms);
        unit=get(Handles.unitsmenu,'Value');
        unitnumber=length(unitmap);
        mode=get(Handles.waveforms,'UserData');
        if isempty(mode)
            mode=0;
        end
        mode=mod(mode+1,3);
        set(Handles.waveforms,'UserData',mode);
        cla;
        axis manual
        
        switch mode 
            case 0 %Waveform plot
                for i=[setdiff(1:unitnumber,unit) unit] %order assures that current unit is drawn last
                    relevant=find(ismember(WaveformInfo.Unit,unitmap{i}));    
                    if length(relevant)>0
                        
                        %       relevant=find(ismember(WaveformInfo.Unit,unitmap{i}));
                        %       temporary_variable = round(linspace(1, length(relevant), 1e2));
                        %       relevant = relevant(temporary_variable);
                        
                        wh=plot([1:size(WaveformInfo.Waveforms,2)]'*ones(1,length(relevant))/FileInfo(1).SamplingRate, ClusterInfo.microVperBit*WaveformInfo.Waveforms(relevant,:)','Color',colors(mod(unitmap{i}+2,256),:));
                    end
                end
            case 1 %Template plot
                if isfield(ClusterInfo,'Centers')
                    for i=1:size(ClusterInfo.Centers,1);
                        plot([1:size(ClusterInfo.Centers,2)]/FileInfo(1).SamplingRate,ClusterInfo.microVperBit*ClusterInfo.Centers(i,:),'Color',colors(mod(ClusterInfo.Units(i)+2,256),:));
                    end
                end      
            case 2 %plot of 10 waveforms from each unit
                for i=1:length(unitmap) %order assures that current unit is drawn last
                    relevant=find(ismember(WaveformInfo.Unit,unitmap{i}));
                    if length(relevant)>0
                        relevant=relevant(round(linspace(1,length(relevant),min(10,length(relevant)))));
                        wh=plot([1:size(WaveformInfo.Waveforms,2)]'*ones(1,length(relevant))/FileInfo(1).SamplingRate, ClusterInfo.microVperBit*WaveformInfo.Waveforms(relevant,:)','Color',colors(mod(unitmap{i}+2,256),:),'LineWidth',2);
                    end
                end
        end
        
        
        
    case 'change_PC'
        set(gcbf,'CurrentAxes',Handles.PCA);
        unit=get(Handles.unitsmenu,'Value');
        unitnumber=length(unitmap);
        mode=get(Handles.PCA,'UserData');
        if isempty(mode)
            mode=0;
        end
        mode=mod(mode+1,3);
        set(Handles.PCA,'UserData',mode);
        switch mode 
            case 0 %Waveform plot
                PCs=[1 2];
            case 1 %Template plot
                PCs=[1 3];     
            case 2 %plot of 10 waveforms from each unit
                PCs=[2 3];
        end
        cla;
        axis auto
        set(gca,'FontSize',15);
        [pc,score,latent] = princomp(WaveformInfo.Waveforms);
        for i=[setdiff(1:unitnumber,unit) unit] %order assures that current unit is drawn last
            relevant=find(ismember(WaveformInfo.Unit,unitmap{i}));    
            if length(relevant)>0
                wh=plot(score(relevant,PCs(1)),score(relevant,PCs(2)), '.','MarkerSize',4, 'Color',colors(mod(unitmap{i}+2,256),:));
            end
        end
        zero=(zeros(size(mean(WaveformInfo.Waveforms)))-mean(WaveformInfo.Waveforms))*pc(:,PCs);
        plot(zero(1),zero(2),'k*','MarkerSize',8)
        xlabel(['PC ' num2str(PCs(1))]);
        ylabel(['PC ' num2str(PCs(2))]);
        if  isfield(ClusterInfo,'Centers')
            clusternum=size(ClusterInfo.Centers,1);
            hold on     
            for i=1:clusternum
                sig1=pc(:,PCs)'*(ClusterInfo.Sigma{i}*pc(:,PCs));%rotate covariance to pc dimensions
                [V,D]=eig(sig1);
                mu=(ClusterInfo.Centers(i,:)-mean(WaveformInfo.Waveforms))*pc(:,PCs);
                ph1(i)=sac_ellipse(2*sqrt(D(1,1)),2*sqrt(D(2,2)),angle(V(1,1)+1i*V(2,1)),mu(1),mu(2));
                set(ph1(i),'LineWidth',2,'Color',colors(mod(ClusterInfo.Units(i)+2,256),:));
            end
        end
        
        
    case 'exit'
        close all; 
end