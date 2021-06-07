function sac_apply_params

global ClusterInfo;

store=ClusterInfo;%save ClusterInfo in temporary variable


%read parameter file
[param_filename,param_pathname]=uigetfile('*ruleset.mat', 'Choose Ruleset file');
if param_filename==0
    return;
end
load ([param_pathname param_filename]);

h=waitbar(0,'Writing Channels...','Position',[250 400 270 60]);

for Channel=1:length(ruleset)% Loop across channels
    waitbar(Channel/length(ruleset),h);
    %if length(setdiff(ruleset(Channel).Units,[0 255]))>0; 
    if length(ruleset(Channel).Units)>0 % assuming you want to reclassify noise only channels
        ClusterInfo=ruleset(Channel);
        ClusterInfo.ChannelNumber=Channel;
        sac_write;
    end
end

ClusterInfo=store;
close(h);