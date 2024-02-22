function[H,L,Length]=same_cut_histogram(sorted_errors)
L=[1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1];%%分割长度的分配 L不能为空矩阵
n1=numel(L);
n=numel(sorted_errors);
Length=zeros(1,n1+1);
for i=1:n1
Length(i+1)=uint32(n*sum(L(1:i))/sum(L));
end
H=zeros(numel(L),511);
for i=1:numel(L)
H(i,:)=hist(sorted_errors(Length(i)+1:Length(i+1)),-255:255);
end