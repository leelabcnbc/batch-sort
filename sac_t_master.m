function unit_num=sac_t_master()

global WaveformInfo
global ClusterInfo
global FileInfo
global Handles
method=['t'];
clear ClusterInfo.Units 

%Initializing

[n,d_wave]=size(WaveformInfo.Waveforms);
options.disp=0;
Penalty_factors=[0.5 0.7 0.9 1 1.2 1.5 2];
switch get(Handles.PCA_num,'Value')
case 1
   d=2;
   N=30;
   PCA=1;
case 2
   d=3;
   N=45;
   PCA=1;
case 3
   d=5;
   N=60;
   PCA=1;
case 4
   d=10;
   N=80;
   PCA=1;
case 5
   d=20;
   N=150;
   PCA=1;
case 6
   d=d_wave;
   N=300;
   PCA=0;
end
N=N/Penalty_factors(get(Handles.Penalty,'Value'));
g=10;

while round (n/N)<g %This seems to make a basic adjustment for the case of a very small # of waveforms
    N=N/2;
    g=g-2;
end

[pc,score] = princomp(WaveformInfo.Waveforms); 
if PCA==1
   global store 
   store=WaveformInfo.Waveforms;
   WaveformInfo.Waveforms=score(:,1:d);
end

size(WaveformInfo.Waveforms)

%Fuzzy c-means sorting
[mu,U]=fcm(WaveformInfo.Waveforms,g,[2,20,1,0]);%20 fuzzy c-means iterations


%estimate error covariance matrix
rep=reshape(repmat(1:g,n,1),g*n,1);
rep_data=repmat(WaveformInfo.Waveforms,g,1);
diffs=rep_data-mu(rep,:);%Distances to cluster center
clear rep_data;
U=U';
for i=1:g
   Sigma{i}=(((U(:,i)*ones(1,d)).*diffs(find(rep==i),:))'*diffs(find(rep==i),:))/sum(U(:,i));
end
%Sigma=cov(WaveformInfo.Waveforms(SortIndices,:))/(g);
Pj=[(sum(U))/sum(sum(U))];

%EM sorting
nu=20;%initial nu value
[Pj,mu,Sigma,nu,Z,U,M]=sac_t_algorithm(Pj,mu,Sigma,nu,N,10,'regular');%10 iterations of regular EM
[Pj,mu,Sigma,nu,Z,U,M]=sac_t_algorithm(Pj,mu,Sigma,nu,N,500,'agglomerate');%modified Mario algorithm
[Pj,mu,Sigma,nu,Z,U,M]=sac_t_algorithm(Pj,mu,Sigma,nu,N,10,'regular');%10 iterations of regular EM
g=length(Pj);
ZU=Z.*U;
Max_M=sac_distribution_inverse(nu,d,0.998);%criterion for allowing 99.8% of spikes
outliers=find(min(M')>Max_M);

if PCA==1 %calculate means and covariances in full waveforms space
   WaveformInfo.Waveforms=store;
   mu=mu*(pc(:,1:d)')+ones(size(mu,1),1)*mean(store);
   d=d_wave;
   for i=1:g
      diffs=WaveformInfo.Waveforms-ones(n,1)*mu(i,:);
      Sigma{i}=(((ZU(:,i)*ones(1,d)).*diffs)'*diffs)/sum(ZU(:,i));;
   end  
end

% [Pj,mu,Sigma,nu,Z,U,M]=sac_t_algorithm(Pj,mu,Sigma,nu,N,1,'regular');%1 iteration of regular EM - full Waveform
% I removed this line since it creates badly conditioned covariance matrices when there is little data.
% How to move from PCA space to full-waveform space needs to be better investigated.
ClusterInfo.Centers=mu;
ClusterInfo.Sigma=Sigma;
ClusterInfo.Proportions=Pj;
[Y,Classification]=max((Z)',[],1);




ClusterInfo.Units=[];
if g>1 %decide which clusters contain local field potentials or overlaps
   for i=1:g
      sig1=pc(:,1:2)'*(ClusterInfo.Sigma{i}*pc(:,1:2));
      f(i)=prod(diag(sig1));
   end
   garbage=find(f>10*median(f)); % Large covariance indicates a garbage collector
   not_garbage=setdiff(1:g,garbage);
   [mu_energy,new_order]=sort(sum(detrend(ClusterInfo.Centers(not_garbage,:)').^2));%sort by ascending energy
   sig1=pc(:,1:2)'*(ClusterInfo.Sigma{not_garbage(new_order(1))}*pc(:,1:2));
   mu1=ClusterInfo.Centers(not_garbage(new_order(1)),:)*pc(:,1:2);
   if sum(mu1.^2)<4*sum(diag(sig1)) % LFP decision
      unit_num=length(new_order)-1;
      ClusterInfo.Units(not_garbage(new_order))=0:unit_num;
   else
      unit_num=length(new_order);
      ClusterInfo.Units(not_garbage(new_order))=1:unit_num;
   end
   ClusterInfo.Units(garbage)=255;
else
   unit_num=1;
   ClusterInfo.Units=1;
end

Classification=ClusterInfo.Units(Classification);

% Classification(find(max((Z.*U)',[],1)<0.5))=255; %outliers are assigned to noise cluster
Classification(outliers)=255;
WaveformInfo.Unit=Classification;

if method==['t']
   ClusterInfo.nu=nu;
else
   ClusterInfo.nu=inf;
end
disp(['nu=' num2str(ClusterInfo.nu)])
