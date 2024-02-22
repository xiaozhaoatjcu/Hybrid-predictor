function[H,L,Length]=per_cut_histogram(sorted_errors,sortsurS)
%%至少三个直方图
%%等分预测误差
feature=sortsurS(end-1,:);%%预测误差的倒数第二行一般是指从小到大排序特征
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
error('直方图数量过大');
end
n=numel(sorted_errors);
Length=zeros(1,numel(L)+1);
H=zeros(numel(L),511);
for i=1:numel(L)%%对直方图进行分段
     Length(i+1)=ceil(n*L(i));
     tmp_H=hist(sorted_errors(Length(i)+1:Length(i+1)),-255:255);
     H(i,:)=tmp_H;
end