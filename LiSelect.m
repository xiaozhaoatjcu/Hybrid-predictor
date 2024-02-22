function[Map_PZ,success,Length,MP,BIGDATA,S_Rate,S_dist]=LiSelect(rate,H,L,Length)
%%作用：
%%安装Li的方法对预测误差分配，并且分配容量
%%安装Li的方法选峰值点P

%%输入：
%%sorted_errors:排序的预测误差
%%rate：嵌入的容量大小
%%输出：

%%Map_PZ输出的峰值点和零值点
%%success嵌入成功的标志
%%Length分割的长度
%%MP根据每段分割的预测误差的峰值点和零值点PZ生出每段预测误差的累计平移量
success=1;
k=8;    %%默认d的范围为0到(k-1)，a=-d-1;
Map_PZ=zeros(numel(L),5);%%默认在预测误差中每个分配片段取2对峰值点和零值点
S_Rate=zeros(1,numel(L));
S_dist=zeros(1,numel(L));
for i=1:numel(L)
[~,tmpzeropos]=find(H(i,:));
Map_PZ(i,3)=tmpzeropos(1)-1-256;%%使用256左边的最右非零数
Map_PZ(i,4)=tmpzeropos(end)+1-256;%%使用256右边的最左非零数
end
j=1;%%j的范围为0-7
m0=0;
tmp_Map_PZ=Map_PZ;
[Map_PZ]=Select(k,L,rate,Map_PZ,H,j,tmp_Map_PZ,m0);
number=size(Map_PZ,2)/5-1;%%多少种可能性
BIGDATA=Map_PZ;
MP=[];
if number==0
success=0;
return
end

D=zeros(1,number);%%存储每种可能性带来的误差
for i=1:number
D(i)=sum(Map_PZ(:,(i+1)*5));
end
[~,pos]=min(D);
Map_PZ=Map_PZ(:,5*pos+1:(pos+1)*5);
for i=1:numel(L)
if Map_PZ(i,5)==0
    Map_PZ(i,1:4)=255;
else
    Map_PZ(i,5)=2;%%存放误差的位置改放峰值点零值点对数2
end
end

for i=1:numel(L)
 if Map_PZ(i,end)==0
 else
     layer=Map_PZ(i,end);
     P=Map_PZ(i,1:layer)+256;
     Z=Map_PZ(i,layer+1:layer+layer)+256;
     [S_Rate(i),S_dist(i)] = rateAndDist(P,Z,H(i,:));
 end
end

for i=1:numel(L)
    
    if Map_PZ(i,5)==0%%在某一段分割预测误差生成的直方图没有用到PZ点，即PZ点的对数为0对
 tmpMP=zeros(511,2);
 
 for j = -255:255
    tmpMP(j+256,:) = tmpMP(j+256,:) + j;                                         
 end

    else
 tmpMP = getModifyMap(Map_PZ(i,1:2),Map_PZ(i,3:4));
    end
 MP=[MP tmpMP];
 
end
%%
function[Map_PZ]=Select(k,L,rate,Map_PZ,H,j,tmp_Map_PZ,m0)
for m1=m0:k-1
tmpj=j+1; 
a=H(j,m1+256);
b=H(j,-m1-1+256);
t_rate=rate-a-b;
tmp_Map_PZ(j,2)=m1;
tmp_Map_PZ(j,1)=-m1-1;%%
[~,D]=rateAndDist(tmp_Map_PZ(j,1:2)+256,tmp_Map_PZ(j,3:4)+256,H(j,:));
tmp_Map_PZ(j,5)=D;
if j<numel(L)
    if t_rate>0
[Map_PZ]=Select(k,L,t_rate,Map_PZ,H,tmpj,tmp_Map_PZ,m1);
    end
end
if t_rate<=0
Map_PZ=[Map_PZ tmp_Map_PZ];
end
end
%%
function [mp] = getModifyMap(P,Z)
%根据每一层平移和嵌入的叠加构造差值映射表
%输入：
%   P，峰值点
%   Z，零点
%输出：
%   mp，差值映射表，若差值i不是峰值点，则mp（i+256，1） = mp（i+256，2）为平移后差值,
%       否则mp（i+256，1）为嵌入0后的差值，mp（i+256，2）为嵌入1后的差值

mp = zeros(511,2);
P = P+256;
Z = Z+256;
for i=1:numel(P)
    if P(i) < Z(i)
        if P(i)+1 <= Z(i)-1
            mp(P(i)+1:Z(i)-1,:) = mp(P(i)+1:Z(i)-1,:) + 1;                 %平移
        end
        mp(P(i),2) = mp(P(i),2) + 1;                                       %嵌入
    else
        if P(i)-1 >= Z(i)+1
            mp(Z(i)+1:P(i)-1,:) = mp(Z(i)+1:P(i)-1,:) - 1;                 %平移
        end
        mp(P(i),2) = mp(P(i),2) - 1;                                       %嵌入
    end
end
for i = -255:255
    mp(i+256,:) = mp(i+256,:) + i;                                         %预测误差加上修改量后的预测误差
end
