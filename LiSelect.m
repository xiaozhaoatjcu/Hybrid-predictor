function[Map_PZ,success,Length,MP,BIGDATA,S_Rate,S_dist]=LiSelect(rate,H,L,Length)
%%���ã�
%%��װLi�ķ�����Ԥ�������䣬���ҷ�������
%%��װLi�ķ���ѡ��ֵ��P

%%���룺
%%sorted_errors:�����Ԥ�����
%%rate��Ƕ���������С
%%�����

%%Map_PZ����ķ�ֵ�����ֵ��
%%successǶ��ɹ��ı�־
%%Length�ָ�ĳ���
%%MP����ÿ�ηָ��Ԥ�����ķ�ֵ�����ֵ��PZ����ÿ��Ԥ�������ۼ�ƽ����
success=1;
k=8;    %%Ĭ��d�ķ�ΧΪ0��(k-1)��a=-d-1;
Map_PZ=zeros(numel(L),5);%%Ĭ����Ԥ�������ÿ������Ƭ��ȡ2�Է�ֵ�����ֵ��
S_Rate=zeros(1,numel(L));
S_dist=zeros(1,numel(L));
for i=1:numel(L)
[~,tmpzeropos]=find(H(i,:));
Map_PZ(i,3)=tmpzeropos(1)-1-256;%%ʹ��256��ߵ����ҷ�����
Map_PZ(i,4)=tmpzeropos(end)+1-256;%%ʹ��256�ұߵ����������
end
j=1;%%j�ķ�ΧΪ0-7
m0=0;
tmp_Map_PZ=Map_PZ;
[Map_PZ]=Select(k,L,rate,Map_PZ,H,j,tmp_Map_PZ,m0);
number=size(Map_PZ,2)/5-1;%%�����ֿ�����
BIGDATA=Map_PZ;
MP=[];
if number==0
success=0;
return
end

D=zeros(1,number);%%�洢ÿ�ֿ����Դ��������
for i=1:number
D(i)=sum(Map_PZ(:,(i+1)*5));
end
[~,pos]=min(D);
Map_PZ=Map_PZ(:,5*pos+1:(pos+1)*5);
for i=1:numel(L)
if Map_PZ(i,5)==0
    Map_PZ(i,1:4)=255;
else
    Map_PZ(i,5)=2;%%�������λ�øķŷ�ֵ����ֵ�����2
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
    
    if Map_PZ(i,5)==0%%��ĳһ�ηָ�Ԥ��������ɵ�ֱ��ͼû���õ�PZ�㣬��PZ��Ķ���Ϊ0��
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
%����ÿһ��ƽ�ƺ�Ƕ��ĵ��ӹ����ֵӳ���
%���룺
%   P����ֵ��
%   Z�����
%�����
%   mp����ֵӳ�������ֵi���Ƿ�ֵ�㣬��mp��i+256��1�� = mp��i+256��2��Ϊƽ�ƺ��ֵ,
%       ����mp��i+256��1��ΪǶ��0��Ĳ�ֵ��mp��i+256��2��ΪǶ��1��Ĳ�ֵ

mp = zeros(511,2);
P = P+256;
Z = Z+256;
for i=1:numel(P)
    if P(i) < Z(i)
        if P(i)+1 <= Z(i)-1
            mp(P(i)+1:Z(i)-1,:) = mp(P(i)+1:Z(i)-1,:) + 1;                 %ƽ��
        end
        mp(P(i),2) = mp(P(i),2) + 1;                                       %Ƕ��
    else
        if P(i)-1 >= Z(i)+1
            mp(Z(i)+1:P(i)-1,:) = mp(Z(i)+1:P(i)-1,:) - 1;                 %ƽ��
        end
        mp(P(i),2) = mp(P(i),2) - 1;                                       %Ƕ��
    end
end
for i = -255:255
    mp(i+256,:) = mp(i+256,:) + i;                                         %Ԥ���������޸������Ԥ�����
end
