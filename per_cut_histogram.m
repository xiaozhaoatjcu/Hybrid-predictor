function[H,L,Length]=per_cut_histogram(sorted_errors,sortsurS)
%%��������ֱ��ͼ
%%�ȷ�Ԥ�����
feature=sortsurS(end-1,:);%%Ԥ�����ĵ����ڶ���һ����ָ��С������������
if feature(1)<0
feature=feature+abs(feature(1));
end
start=feature(end);
key=0;
L=[];
while ~key
TMP=sum(feature<=start);
start=start/1.28;
L=[TMP L];
    if  TMP<8000
        key=1;
    end
end
L=L./numel(sorted_errors);
if numel(L)>31
error('ֱ��ͼ��������');
end
n=numel(sorted_errors);
Length=zeros(1,numel(L)+1);
H=zeros(numel(L),511);
for i=1:numel(L)%%��ֱ��ͼ���зֶ�
     Length(i+1)=ceil(n*L(i));
     tmp_H=hist(sorted_errors(Length(i)+1:Length(i+1)),-255:255);
     H(i,:)=tmp_H;
end