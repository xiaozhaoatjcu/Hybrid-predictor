function [MP,mn,P,Z] = newPickPeak(H,layer,rate)
%枚举选择零极点并计算失真,点对数已指定
%输入：
%   H，直方图
%   layer，点对数
%   rate，嵌入容量
%输出：
%   mn，最小失真
%   P，峰值点
%   Z，零点

%对直方图频率由大到小排序，选取前MaxNump个作为备选的peak点集
[tH, tpos] = sort(H,'descend'); %排序频率值由高到低
numz = sum(tH == 0);%零值点的对数
nump = sum(tH >18);%峰值点的对数
MaxNump = layer;
for i = ceil(layer/2)*2:2:ceil(nump/2)*2   %如果nump/2>layer/2得到的将是一个空矩阵
    if factorial(i)/factorial(layer)/factorial(i-layer) < 10^6   %C(m，n)<10^6
        MaxNump = i;
    end
end
nump = MaxNump;

pSet=tpos(1:nump);%初始化峰值零点集合
zSet=sort(tpos(numel(H)-numz+1:numel(H)));
pset = tH(1:nump);

%计算所有满足容量且刚好为layer个peak点的组合
MP = pick(rate,layer,pset);
[row,col] = size(MP);

if row ==0
    mn = 0;
    P =[];
    Z = [];
    return;
end

%对每个peak点组合，为其选择zero点组合并计算失真
Pst = sum(zSet<256)+1;
D = zeros(1,row);
for i = 1:row
    [D(i),P,Z] = selectZAndCalDist(H,MP(i,:),tpos,zSet,Pst,nump);
end
%选择失真最小的一组并选择zero点
[mn,idx] = min(D);
[mn,P,Z] = selectZAndCalDist(H,MP(idx,:),tpos,zSet,Pst,nump);
P = P-256;
Z = Z-256;
%%
function [M] = pick(rate,layer,pSet)
%对给定的peak点集合，层数以及容量，计算可能的peak点组合
%输入：
%   rate，容量
%   layer，点对数
%   pSet，峰值点频率集合
%输出：
%   M，可行的峰值点组合

%numLeft为pSet大小
numLeft = numel(pSet);
%若pSet大小小于层数，则无可行peak点组合，返回
if numLeft < layer 
    M = [];
    return;
end
%若pSet频率和小于容量，则无可行peak点组合，返回
if sum(pSet(1,1:layer)) < rate
    M = [];
    return;
end

if layer ==1
    %若只剩最后一层，则在pSet中选择所有count个满足的频率点，返回count种组合
    count = numLeft;
    for i = 1:numLeft
        if pSet(1,i) < rate
            count = i-1;
            break;
        end
    end
    if count == 0
        M = [];
        return;
    elseif count == numLeft
        M = eye(count);
    else
        M = [eye(count),zeros(count,numLeft-count)];
    end
else
    %若剩下不止一层，则搜索小于容量rate的频率点
    startP = 1;
    while startP <= numLeft && pSet(startP) >= rate
        startP  = startP + 1;%剔除掉pSet中所有单个高度大于rate的点
    end
    %若小于rate的频率点数小于layer，则无可行peak点组合，返回
    if startP > numLeft || numLeft - startP + 1 < layer
        M = [];
        return;
    else
        %否则对第一个小于rate的频率点，选择或者不选择此点，然后对此函数递归
        M1 = pick(rate,layer,pSet(1,startP+1:numLeft));%不选择此点
        M2 = pick(rate-pSet(1,startP),layer-1,pSet(1,startP+1:numLeft));%选择此点
        %返回所有组合M
        [r1,c1] = size(M1);
        [r2,c2] = size(M2);
        if r1 == 0
            if r2 == 0
                M = [];
            else
                M = [zeros(r2,startP-1),ones(r2,1),M2];
            end
        else
            if r2 == 0
                M = [zeros(r1,startP),M1];
            else
                M = [zeros(r1,startP),M1;zeros(r2,startP-1),ones(r2,1),M2];
            end
        end
    end
end
M = boolean(M);%布尔数学变量
%%
function [dist,P,Z] = selectZAndCalDist(H,FP,tpos,zSet,Pst,nump)
%对每个峰值点组合，选择对应的零点并计算失真
%输入：
%   H，直方图
%   FP，峰值点使用标志，FP（i）为1表示使用频率第i大的峰值点
%   tpos，tpos（i）表示第i大频率的峰值点
%   zSet，零点集合
%   Pst，正的零点在zSet中的开始索引
%   nump，最大点对数
%输出：
%   dist，失真
%   P，峰值点
%   Z，零点
L = sum(FP>0);
P = zeros(1,L);
Z = zeros(1,L);
k = 1;
for j = 1:nump
    if FP(j) == 1
        P(k) = tpos(j);
        k = k+1;
    end
end

idP = Pst;
idN= Pst - 1;
for i=1:L
    if(P(i)>255)
        Z(i) = zSet(idP);
        idP = idP+1;
    else
        Z(i) = zSet(idN);
        idN = idN-1;
    end
end
[rate,dist] = rateAndDist(P,Z,H);