function[H_definite_infor]=simulate_entropy(Length,sorted_errors)
m=numel(Length)-1;
High=zeros(1,m);
Wide=zeros(1,m);
entr=zeros(1,m); 
H=zeros(m,511);
for i=1:m
tmp_H=sorted_errors((Length(i)+1):Length(i+1));
H(i,:)=hist(tmp_H,-255:255);
tmpdata=find(H(i,:));
Wide(i)=tmpdata(end)-tmpdata(1);
High(i)=max(H(i,:));
entr(i)=entropy(sorted_errors((Length(i)+1):Length(i+1)));
end
H_definite_infor.High=High;
H_definite_infor.entr=entr;
H_definite_infor.Wide=Wide;
H_definite_infor.H=H;