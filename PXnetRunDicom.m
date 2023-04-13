clear all

load PXnet.mat
load('Case.mat')
%% import your CT
% folder=uigetdir('../data','Select directory containing CT images');
% [CTcont, VS] = dicomreadfolder2(folder);
% 
% if VS(3)~=3
%     CTcontB=imresize3(double(CTcont),[size(CTcont,1) size(CTcont,2) round(size(CTcont,3)/(3/VS(3)))]);
% else
%     CTcontB=double(CTcont);
% end
%%
CTcontB=double(CTcontB);
if size(CTcontB,3)<48
    CTcontB(:,:,end+1:48)=0;
end
    
%0-2048
CT=((CTcontB+(VS(4)+1024)))-1024;
thr1=0;
thr2=490;
CTthr=CT;
c1=find(CT<thr1+0);
c2=find(CT>thr2+0);

CTthr(c1)=0;
CTthr(c2)=0;


N=size(CTthr,3);
Nres=mod(N,48);
Nadd=48-Nres;
CTthrAdd=cat(3,CTthr,zeros(512,512,Nadd));

N48=size(CTthrAdd,3)/48;

for n=1:N48
    CTcrop=CTthrAdd(:,:,48*(n-1)+1:48*n);
    C = semanticseg(CTcrop,net);
    yerM=find(C=='M');
    yerC=find(C=='C');
    Cout=zeros(512,512,48);
    Cout(yerC)=1;
    Cout(yerM)=2;
    
    CTout1(:,:,48*(n-1)+1:48*n)=Cout;
end

for n=1:N48-1
    CTcrop=CTthrAdd(:,:,48*(n-1)+1+24:48*n+24);
    C = semanticseg(CTcrop,net);
    yerM=find(C=='M');
    yerC=find(C=='C');
    Cout2=zeros(512,512,48);
    Cout2(yerC)=1;
    Cout2(yerM)=2;
    
    CTout2(:,:,48*(n-1)+1+24:48*n+24)=Cout2;
end

CTout1(:,:,end-Nadd+1:end)=[];
diflngth=size(CTout1,3)-size(CTout2,3);
if diflngth>0
    CTout2=cat(3,CTout2,zeros(512,512,diflngth));
else
    CTout2(:,:,end+diflngth+1:end)=[];
end


%% combine windos
CTout=zeros(512,512,size(CTout1,3));
indexC=find(CTout1==1 | CTout2==1); 
CTout(indexC)=1;
indexM=find(CTout1==2 | CTout2==2);
CTout(indexM)=2;

%% selecting kidney
CTkn=CTout;

CTc=zeros(512,512,size(CTkn,3))>0;
CTm=zeros(512,512,size(CTkn,3))>0;

mc=find(CTkn==1);
mm=find(CTkn==2);
CTc(mc)=1;
CTm(mm)=1;
se1=strel('sphere',1);
se3=strel('sphere',3);
se10=strel('disk',10);

CTc1=imopen(CTc,se1);
CTm1=imopen(CTm,se1);

CTc1=imclose(CTc1,se10);
CTc1=imfill(CTc1,8,'holes');
CTm1=imdilate(CTm1,se10);
CTm1=imdilate(CTm1,se3);

CTcm=CTc1+CTm1;
CTcm1=CTcm==2;
LL=bwlabeln(CTcm1,26);

stats=regionprops('table',LL,'area','centroid');
aa=table2array(stats(:,1));
[areaSort,yer]=sort(aa,'descend');
cc=table2array(stats(yer(1),2));
center(1,:)=cc;
cc=table2array(stats(yer(2),2));
center(2,:)=cc;
LM1=LL==yer(1);
LM2=LL==yer(2);
CTMedulla=CTkn==2;
LM=CTMedulla & (LM1 | LM2);


% centroids of PXNetnet medullas
CTCortex=CTkn==1;
L=bwlabeln(CTCortex,26);
stats=regionprops('table',L,'area','centroid');
aa=table2array(stats(:,1))*VS(1)*VS(2)*3/1000; %ml
buyukCortex=find(aa>50);
center2=table2array(stats(:,2));
center21=center2(buyukCortex,:);

%Comperison of Euclidian Distances
for k=1:2
    for m=1:size(center21,1)
        dst(m,k)=sum((center(k,:)-center21(m,:)).^2);
    end
end

[~,yer1]=min(dst(:,1));
[~,yer2]=min(dst(:,2));
LC1=L==buyukCortex(yer1);
LC2=L==buyukCortex(yer2);

LCM=uint8(LC1+LC2+LM*2);





