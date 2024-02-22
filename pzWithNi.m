function [P,Z] = pzWithNi(H,rate)
%Ni ZhiCheng���˷�����ѡ��
%���룺
%   H��ֱ��ͼ
%   rate������
%�����
%   P����ֵ�㣬-255~255��Χ��
%   Z����㣬-255~255��Χ��
[fsort,idx] = sort(H,'descend');
lv = 0;
sum = 0;
while sum < rate
    lv = lv + 1;
    sum = sum + fsort(lv);
    if lv > 40
        error('Rate is too large');
    end
end
P = idx(1,1:lv);
Z = zeros(1,lv);
for i = 1:lv
    offset = 1;
    while H(P(i) - offset) && H(P(i) + offset)
        offset = offset + 1;
    end
    if H(P(i) - offset)
        Z(i) = P(i) + offset;
    else
        Z(i) = P(i) - offset;
    end
    H(Z(i)) = 1;
end
P = P - 256;
Z = Z - 256;