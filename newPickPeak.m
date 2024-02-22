function [MP,mn,P,Z] = newPickPeak(H,layer,rate)
%ö��ѡ���㼫�㲢����ʧ��,�������ָ��
%���룺
%   H��ֱ��ͼ
%   layer�������
%   rate��Ƕ������
%�����
%   mn����Сʧ��
%   P����ֵ��
%   Z�����

%��ֱ��ͼƵ���ɴ�С����ѡȡǰMaxNump����Ϊ��ѡ��peak�㼯
[tH, tpos] = sort(H,'descend'); %����Ƶ��ֵ�ɸߵ���
numz = sum(tH == 0);%��ֵ��Ķ���
nump = sum(tH >18);%��ֵ��Ķ���
MaxNump = layer;
for i = ceil(layer/2)*2:2:ceil(nump/2)*2   %���nump/2>layer/2�õ��Ľ���һ���վ���
    if factorial(i)/factorial(layer)/factorial(i-layer) < 10^6   %C(m��n)<10^6
        MaxNump = i;
    end
end
nump = MaxNump;

pSet=tpos(1:nump);%��ʼ����ֵ��㼯��
zSet=sort(tpos(numel(H)-numz+1:numel(H)));
pset = tH(1:nump);

%�����������������Ҹպ�Ϊlayer��peak������
MP = pick(rate,layer,pset);
[row,col] = size(MP);

if row ==0
    mn = 0;
    P =[];
    Z = [];
    return;
end

%��ÿ��peak����ϣ�Ϊ��ѡ��zero����ϲ�����ʧ��
Pst = sum(zSet<256)+1;
D = zeros(1,row);
for i = 1:row
    [D(i),P,Z] = selectZAndCalDist(H,MP(i,:),tpos,zSet,Pst,nump);
end
%ѡ��ʧ����С��һ�鲢ѡ��zero��
[mn,idx] = min(D);
[mn,P,Z] = selectZAndCalDist(H,MP(idx,:),tpos,zSet,Pst,nump);
P = P-256;
Z = Z-256;
%%
function [M] = pick(rate,layer,pSet)
%�Ը�����peak�㼯�ϣ������Լ�������������ܵ�peak�����
%���룺
%   rate������
%   layer�������
%   pSet����ֵ��Ƶ�ʼ���
%�����
%   M�����еķ�ֵ�����

%numLeftΪpSet��С
numLeft = numel(pSet);
%��pSet��СС�ڲ��������޿���peak����ϣ�����
if numLeft < layer 
    M = [];
    return;
end
%��pSetƵ�ʺ�С�����������޿���peak����ϣ�����
if sum(pSet(1,1:layer)) < rate
    M = [];
    return;
end

if layer ==1
    %��ֻʣ���һ�㣬����pSet��ѡ������count�������Ƶ�ʵ㣬����count�����
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
    %��ʣ�²�ֹһ�㣬������С������rate��Ƶ�ʵ�
    startP = 1;
    while startP <= numLeft && pSet(startP) >= rate
        startP  = startP + 1;%�޳���pSet�����е����߶ȴ���rate�ĵ�
    end
    %��С��rate��Ƶ�ʵ���С��layer�����޿���peak����ϣ�����
    if startP > numLeft || numLeft - startP + 1 < layer
        M = [];
        return;
    else
        %����Ե�һ��С��rate��Ƶ�ʵ㣬ѡ����߲�ѡ��˵㣬Ȼ��Դ˺����ݹ�
        M1 = pick(rate,layer,pSet(1,startP+1:numLeft));%��ѡ��˵�
        M2 = pick(rate-pSet(1,startP),layer-1,pSet(1,startP+1:numLeft));%ѡ��˵�
        %�����������M
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
M = boolean(M);%������ѧ����
%%
function [dist,P,Z] = selectZAndCalDist(H,FP,tpos,zSet,Pst,nump)
%��ÿ����ֵ����ϣ�ѡ���Ӧ����㲢����ʧ��
%���룺
%   H��ֱ��ͼ
%   FP����ֵ��ʹ�ñ�־��FP��i��Ϊ1��ʾʹ��Ƶ�ʵ�i��ķ�ֵ��
%   tpos��tpos��i����ʾ��i��Ƶ�ʵķ�ֵ��
%   zSet����㼯��
%   Pst�����������zSet�еĿ�ʼ����
%   nump���������
%�����
%   dist��ʧ��
%   P����ֵ��
%   Z�����
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