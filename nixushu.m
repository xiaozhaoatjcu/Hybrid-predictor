function[IN,nixushutime]=nixushu(err)
%%ͨ����״����������
%%err:��������
%%IN:������

tt1=clock();
m=numel(err);
temp=zeros(1,m);
[a,tpos]=sort(err,'descend');
IN=0;
for i=1:m
    temp(1,tpos(i))=1;
    IN=IN+sum(temp(1,1:tpos(i)-1));
    i=i+1;
end
i=1;
same=0;
extrosame=1;
for i=1:m-1    %%��ȥ��ͬ����Ŀ����
    if a(i)==a(i+1);
        same=same+extrosame;
        extrosame=extrosame+1;
    else
    extrosame=1;
    end
    i=i+1;
end
    IN=IN-same;
 tt2=clock();
 nixushutime=etime(tt2,tt1);