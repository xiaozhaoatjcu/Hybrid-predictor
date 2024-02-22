function [length,layer,P,Z,idx,Dist,DS,minPer,calCount] = devideAndPick2(err,rate,maxPer)
%搜索最优分割位置,采用类似二分搜索
%输入：
%   err，排序后的差值
%   rate，嵌入容量
%   maxPer，限定的最大分割百分比
%输出：
%   length,分割长度
%   layer，嵌入层数
%   P，峰值点
%   Z，零点
%   idx，最优分割百分比
%   Dist，最优分割点失真
%   DS，所有分割点的失真
%   minPer，最小分割百分比
%   calCount，总共测试的百分比个数

DS = ones(1,100).*100000000;
Layer = ones(1,100);
total = numel(err);
P=[];
Z=[];
Dist = 100000000;
maxLayer = 25;

%对于1%——100%的分割点，生成直方图并进行选点，假定最大的嵌入层数不超过maxLayer层
minPer = ceil(double(rate)/numel(err)*100);%小于minPer的百分比是无法满足嵌入容量的需求的

curArea = ones(4,2);                                                     %当前所有待搜索区域
tmpArea = ones(4,2);                                                     %可能的待搜索区域
calCount = 2;                                                              
curC = 1;                                                                  %当前待搜索区域个数
tmpC = 0;                                                                  %可能的待搜索区域个数

H = hist(err(1,1:uint32(total*minPer*0.01)),-255:255);                     %得到最小分割百分比的直方图；计算最小分割百分比处的失真
[tP,tZ,Layer(1,minPer),DS(1,minPer)] = performOfOneHist(H,rate,maxLayer);
if DS(1,minPer) ~= 0 &&  DS(1,minPer) < Dist
    Dist = DS(1,minPer);
    P = tP;
    Z = tZ;
end
H = hist(err(1,1:uint32(total*maxPer*0.01)),-255:255);                     %计算最大分割百分比处的失真
[tP,tZ,Layer(1,maxPer),DS(1,maxPer)] = performOfOneHist(H,rate,maxLayer);
if DS(1,maxPer) ~= 0 &&  DS(1,maxPer) < Dist
    Dist = DS(1,maxPer);
    P = tP;
    Z = tZ;
end
curArea(1,:) =[minPer,maxPer];                                             %初始的待搜索区域组

if maxPer == minPer+1
    curC = 0;
end

while curC > 0                                                             %待搜索区域不为空时继续
    tmpC = 0;
    for i=1:curC                                                           %对每一个待搜索区域，计算其中间点的失真
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
    for i = 1:tmpC                                                         %从可能的待搜索区域中选择实际需要搜索的区域，替换原待搜索区域
        if tmpArea(i,1) == midx || tmpArea(i,2) == midx
            curC = curC + 1;
            curArea(curC,:) = tmpArea(i,:);
        end
    end
end

%找出最小失真的分割位置
[mn,idx]=min(DS);
length = uint32(total*idx*0.01);
layer = Layer(1,idx);
if mn == 100000000
    length = 1;
end
%%
function [P,Z,layer,Dist] = performOfOneHist(H,rate,maxLayer)
%峰值点-零点选择并计算失真，点对数小时用枚举搜索，否则用GA
%输入：
%   H，直方图
%   rate，容量限制
%   maxLayer，最大点对数
%输出：
%   P，峰值点
%   Z，零点
%   layer，点对数
%   Dist，失真

P = [];
Z = [];
layer = 1;
Dist = 100000000;
[tH,tpos] = sort(H,'descend');
tSum = 0;
minL=0;%最小点对数
for j=1:maxLayer
    tSum = tSum+tH(j);
    if tSum >= rate
        minL = j;
        break;
    end
end%由嵌入容量确定最小点对数
%该直方图容量不足以嵌入时跳到下一个分割点，否则继续
if tSum < rate
    return;
end
maxL = maxLayer;%最大点对数
if maxL > sum(H >18)
    maxL = sum(H>18);
end%最大点对数取sum（H>18）和maxlayer中的较小数
%层数j从最小层数minL开始，计算j层嵌入的最小失真
founded = 0;
for j = minL:7
    [MP,tDist,tP,tZ] = newPickPeak(H,j,rate);%%不再考虑峰值点和零值点本身对嵌入容量的影响
    if tDist == 0
        continue;
    elseif tDist < Dist
        Dist = tDist;
        layer = j;
        P = tP;
        Z = tZ;
    else
        founded = 1;                                                       %若点对数增大失真随之增大，则跳出
        break;
    end
end
