function [P,Z] =  pzWithXuan(H,rate)
%Xuan�ķ�ֵ������Ż�ѡ��
%���룺
%   H,ֱ��ͼ
%   rate,����
%�����
%   P����ֵ�㣬-255~255��Χ��
%   Z����㣬-255~255��Χ��
pp=zeros(1,sum(H(1,256:511)>0));
np=zeros(1,sum(H(1,1:255)>0));
pz=zeros(1,256-numel(pp));
nz=zeros(1,255-numel(np));
maxT=0;
minT=0;
cp = 1;
cz = 1;
for i=255:-1:1
    if H(i) 
        np(cp) = i;
        cp = cp + 1;
    else
        nz(cz) = i;
        cz = cz + 1;
    end
end
cp = 1;
cz = 1;
for i=256:511
    if H(i)
        pp(cp) = i;
        cp = cp + 1;
    else
        pz(cz) = i;
        cz = cz + 1;
    end
end
if 256-nz(1)<pz(1)-256
    maxT=256-nz(1)-1;%��ֵ�����ֵ��֮�����������ѡȡ���н�С��
else
    maxT=pz(1)-256-1;
end
s=H(256);
for i=1:maxT+1
    if s>=rate
        break;
    end
    s=s+H(256+i)+H(256-i);
    minT=i;
end
if minT==maxT+1
    P=0;
    Z=P;
    return
end
tmp = 1;
while H(tmp) == 0
    tmp = tmp+1;
end
tmp = tmp -1;
nz = tmp:-1:1;
tmp = 511;
while H(tmp) == 0;
    tmp = tmp-1;
end
tmp = tmp+1;
pz=tmp:511;
tmd=100000000;
for t=minT:1:maxT
    R=0;
    D=0;
    i=t;
    tP=[];
    tZ=[];
    while R<rate
        R=R+H(i+256);
        if i>=0
            tP=[tP,pp(i+1)];
            if t-i+1 > numel(pz)%i��СΪ0
                P = [254,-254];
                Z = [255,-255];
                return;
            end
            tZ=[tZ,pz(t-i+1)];
            i=-i;
        else
            tP=[tP,np(-i)];
            if t+i+1 > numel(nz)%i��СΪ-1
                P = [254,-254];
                Z = [255,-255];
                return;
            end
            tZ=[tZ,nz(t+i+1)];
            i=-i-1;
        end
    end
    [R,D]=rateAndDist(tP,tZ,H);
    if D<tmd
        tmd=D;
        P=tP;
        Z=tZ;
    end
end
P=P-256;
Z=Z-256;