%% Non-iterative optimal blind deconvolution and its application to machine condition monitoring
% 非迭代最优盲解卷积代码（Code）Code by Liu He
% 期刊: Mechanical Systems and Signal Processing
% 论文下载: doi: 10.1016/j.ymssp.2025.113715
%          https://www.sciencedirect.com/science/article/pii/S0888327025014165?via%3Dihub
% email: aresmiki@163.com, 
% 2025.9  
%%
clc
close all
clear all
%% 加载测试数据集
fileFolder=fullfile('/Users/arsmiki/2nd_test'); %请下载IMS数据集，并替换成自己的目录
dirOutput=dir(fullfile(fileFolder,'*.39'));
fileNames={dirOutput.name};
%% 加载噪声数据（基础振动）
norm_data=load('/Users/arsmiki/x1norm_5k_data.txt');%请下载该噪声数据集，并替换成自己的目录
N=20000;  
norm_data=0.5*norm_data(1:N);
%% 数据截断以及噪声叠加
sN=10000;  %数据截断长度
sL=5000;   %截断数据的重叠长度
ordata=zeros(N,984);
pjdordata=[]; %用于存储叠加噪声后的截断数据集
for i=1:length(fileNames)
    i
    rbt=load(fullfile(fileFolder,fileNames{i}));
    rbt1=rbt(1:N,1)+norm_data;
    ordata(:,i)=rbt1./norm(rbt1,1);
    for j=0:1:floor((N-sN)/sL)
        pjdordata=[pjdordata,rbt1((1+j*sL):(j*sL+sN))];
    end
end
%% 主程序部分
close all
selc=0; %选择0使用DCT，选择1使用FFT
zp=1; % 正负样本大小，选择1意味着正样本1个、负样本1个
bes=size(pjdordata,2)/984; %批次大小，将1个样本通过截断操作后产生的样本大小
np=20; %滤波器长度大小
stat=tic;
%% FFT
if selc==1 
Xh=get_hankel(mean(pjdordata(:,1:zp*bes),2),np);
F=fft(eye(sN-np+1,sN-np+1))*2/(sN-np+1);
F=F(1:floor((sN-np+1)/2),:);
F(1:5,:)=0;
F=[real(F);imag(F)];
rh=[];
obl=[];
for epc=zp+1:984-zp+1
    epc
    Xf=get_hankel(mean(pjdordata(:,(epc-1)*bes+1:(epc+zp-1)*bes),2),np);
    DX=F*(Xf-Xh);
    BX=F*Xh;
    Rxwx=DX'*DX; %沿着时间轴积分
    Rxbx=BX'*BX;
    [h,lambda]=eigs(inv(Rxbx)*Rxwx,1);
    rh=[rh,h];
    % lx = filter(h,1,pjdordata(:,(epc+zp-1)*bes));
    lx = filter(h,1,ordata(:,epc+zp-1));
    bl=abs(fft(abs(hilbert(lx))))*2/length(lx);
    obl=[obl,bl(1:1000)];
end
%% 使用DCT
else
Xh=get_hankel(mean(pjdordata(:,1:zp*bes),2),np);
F=dctmtx(sN-np+1);
F(1:5,:)=0;
rh=[];
obl=[];
for epc=zp+1:984-zp+1
    epc
    Xf=get_hankel(mean(pjdordata(:,(epc-1)*bes+1:(epc+zp-1)*bes),2),np);
    DX=F*(Xf-Xh);
    BX=F*Xh;
    Rxwx=DX'*DX; %沿着时间轴积分
    Rxbx=BX'*BX;
    [h,lambda]=eigs(inv(Rxbx)*Rxwx,1);
    rh=[rh,h];
    % lx = filter(h,1,pjdordata(:,(epc+zp-1)*bes));
    lx = filter(h,1,ordata(:,epc+zp-1));
    bl=abs(fft(abs(hilbert(lx))))*2/length(lx);
    obl=[obl,bl(1:1000)];
end
end
ctime=toc(stat);
%% 归一化
oblp=[];
oblp1=[];
for i=1:size(obl,2)
    i
    tmmq=obl(:,i)./norm(obl(2:end,i),2);
    tmm=obl(:,i)./max(obl(2:end,i));
    oblp=[oblp,tmm];
    oblp1=[oblp1,tmmq];
end
%% 3D绘图
alpha1=20000/length(lx)*(0:1:length(lx)-1);
figure
mesh([1:size(oblp1,2)]+2*zp-1,alpha1(2:1000),oblp1(2:1000,:))
colormap(parula)
clim([0.05 0.2])
view(-40,69)







