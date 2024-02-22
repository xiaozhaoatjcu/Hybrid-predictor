function [length,layer,P,Z,idx,Dist,DS,minPer,calCount] = devideAndPick2(err,rate,maxPer)
%�������ŷָ�λ��,�������ƶ�������
%���룺
%   err�������Ĳ�ֵ
%   rate��Ƕ������
%   maxPer���޶������ָ�ٷֱ�
%�����
%   length,�ָ��
%   layer��Ƕ�����
%   P����ֵ��
%   Z�����
%   idx�����ŷָ�ٷֱ�
%   Dist�����ŷָ��ʧ��
%   DS�����зָ���ʧ��
%   minPer����С�ָ�ٷֱ�
%   calCount���ܹ����Եİٷֱȸ���

DS = ones(1,100).*100000000;
Layer = ones(1,100);
total = numel(err);
P=[];
Z=[];
Dist = 100000000;
maxLayer = 25;

%����1%����100%�ķָ�㣬����ֱ��ͼ������ѡ�㣬�ٶ�����Ƕ�����������maxLayer��
minPer = ceil(double(rate)/numel(err)*100);%С��minPer�İٷֱ����޷�����Ƕ�������������

curArea = ones(4,2);                                                     %��ǰ���д���������
tmpArea = ones(4,2);                                                     %���ܵĴ���������
calCount = 2;                                                              
curC = 1;                                                                  %��ǰ�������������
tmpC = 0;                                                                  %���ܵĴ������������

H = hist(err(1,1:uint32(total*minPer*0.01)),-255:255);                     %�õ���С�ָ�ٷֱȵ�ֱ��ͼ��������С�ָ�ٷֱȴ���ʧ��
[tP,tZ,Layer(1,minPer),DS(1,minPer)] = performOfOneHist(H,rate,maxLayer);
if DS(1,minPer) ~= 0 &&  DS(1,minPer) < Dist
    Dist = DS(1,minPer);
    P = tP;
    Z = tZ;
end
H = hist(err(1,1:uint32(total*maxPer*0.01)),-255:255);                     %�������ָ�ٷֱȴ���ʧ��
[tP,tZ,Layer(1,maxPer),DS(1,maxPer)] = performOfOneHist(H,rate,maxLayer);
if DS(1,maxPer) ~= 0 &&  DS(1,maxPer) < Dist
    Dist = DS(1,maxPer);
    P = tP;
    Z = tZ;
end
curArea(1,:) =[minPer,maxPer];                                             %��ʼ�Ĵ�����������

if maxPer == minPer+1
    curC = 0;
end

while curC > 0                                                             %����������Ϊ��ʱ����
    tmpC = 0;
    for i=1:curC                                                           %��ÿһ�����������򣬼������м���ʧ��
        tmpIdx = floor((curArea(i,1)+curArea(i,2))/2);
        H = hist(err(1,1:uint32(total*tmpIdx*0.01)),-255:255);
        [tP,tZ,Layer(1,tmpIdx),DS(1,tmpIdx)] = performOfOneHist(H,rate,maxLayer);
        if DS(1,tmpIdx) ~= 0 &&  DS(1,tmpIdx) < Dist
            Dist = DS(1,tmpIdx);
            P = tP;
            Z = tZ;
        end
        calCount  = calCount + 1;
        if tmpIdx > curArea(i,1)+1
            tmpC = tmpC+1;
            tmpArea(tmpC,:) = [curArea(i,1),tmpIdx];
        end
        if tmpIdx + 1 < curArea(i,2)
            tmpC = tmpC + 1;
            tmpArea(tmpC,:) = [tmpIdx,curArea(i,2)];
        end 
    end
    [mdist,midx] = min(DS);
    curC = 0;
    for i = 1:tmpC                                                         %�ӿ��ܵĴ�����������ѡ��ʵ����Ҫ�����������滻ԭ����������
        if tmpArea(i,1) == midx || tmpArea(i,2) == midx
            curC = curC + 1;
            curArea(curC,:) = tmpArea(i,:);
        end
    end
end

%�ҳ���Сʧ��ķָ�λ��
[mn,idx]=min(DS);
length = uint32(total*idx*0.01);
layer = Layer(1,idx);
if mn == 100000000
    length = 1;
end
%%
function [P,Z,layer,Dist] = performOfOneHist(H,rate,maxLayer)
%��ֵ��-���ѡ�񲢼���ʧ�棬�����Сʱ��ö��������������GA
%���룺
%   H��ֱ��ͼ
%   rate����������
%   maxLayer���������
%�����
%   P����ֵ��
%   Z�����
%   layer�������
%   Dist��ʧ��

P = [];
Z = [];
layer = 1;
Dist = 100000000;
[tH,tpos] = sort(H,'descend');
tSum = 0;
minL=0;%��С�����
for j=1:maxLayer
    tSum = tSum+tH(j);
    if tSum >= rate
        minL = j;
        break;
    end
end%��Ƕ������ȷ����С�����
%��ֱ��ͼ����������Ƕ��ʱ������һ���ָ�㣬�������
if tSum < rate
    return;
end
maxL = maxLayer;%�������
if maxL > sum(H >18)
    maxL = sum(H>18);
end%�������ȡsum��H>18����maxlayer�еĽ�С��
%����j����С����minL��ʼ������j��Ƕ�����Сʧ��
founded = 0;
for j = minL:7
    [MP,tDist,tP,tZ] = newPickPeak(H,j,rate);%%���ٿ��Ƿ�ֵ�����ֵ�㱾���Ƕ��������Ӱ��
    if tDist == 0
        continue;
    elseif tDist < Dist
        Dist = tDist;
        layer = j;
        P = tP;
        Z = tZ;
    else
        founded = 1;                                                       %�����������ʧ����֮����������
        break;
    end
end
