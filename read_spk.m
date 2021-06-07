fid=fopen('t1.spk.1','r');
s=fread(fid,[4 inf],'int16');
a=reshape(mean(s,1),[80 962480/80]);
global WaveformInfo;
global ClusterInfo;
WaveformInfo.Unit=[];
WaveformInfo.Waveforms=a';
ClusterInfo.microVperBit=1;
