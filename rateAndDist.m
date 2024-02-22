function [rate,dist] = rateAndDist(P,Z,H)
% 计算多点（非多层）嵌入的嵌入容量以及失真
%输出：
%   dist:失真
%   rate:容量
%输入：
%   P峰值点（均未漂移的）（在1~511范围内，相当于实际值加256后转换为数组H的索引）
%   Z零值点（均未漂移的）（在1~511范围内，相当于实际值加256后转换为数组H的索引）
%   H原直方图，H(i)代表预测差值i-256的频率
m=numel(P);
n=numel(H);
C=double(zeros(1,n));
PK=double(zeros(1,n));
sig = sign(Z-P);
PK(P) = sig;%给H中的P，Z点位值左右的位置定义正负值，P在Z左边为1，p在Z右边为-1
for i=1:m
    C(1,P(i)+sig(i):sig(i):Z(i)-sig(i))=C(1,P(i)+sig(i):sig(i):Z(i)-sig(i)) + sig(i);%给涉及到平移的每一个点+1,-1
end
rate = sum(H(P));
dist = 0.5*double(sum(H.*(C.^2+(C+PK).^2))); 