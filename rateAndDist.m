function [rate,dist] = rateAndDist(P,Z,H)
% �����㣨�Ƕ�㣩Ƕ���Ƕ�������Լ�ʧ��
%�����
%   dist:ʧ��
%   rate:����
%���룺
%   P��ֵ�㣨��δƯ�Ƶģ�����1~511��Χ�ڣ��൱��ʵ��ֵ��256��ת��Ϊ����H��������
%   Z��ֵ�㣨��δƯ�Ƶģ�����1~511��Χ�ڣ��൱��ʵ��ֵ��256��ת��Ϊ����H��������
%   Hԭֱ��ͼ��H(i)����Ԥ���ֵi-256��Ƶ��
m=numel(P);
n=numel(H);
C=double(zeros(1,n));
PK=double(zeros(1,n));
sig = sign(Z-P);
PK(P) = sig;%��H�е�P��Z��λֵ���ҵ�λ�ö�������ֵ��P��Z���Ϊ1��p��Z�ұ�Ϊ-1
for i=1:m
    C(1,P(i)+sig(i):sig(i):Z(i)-sig(i))=C(1,P(i)+sig(i):sig(i):Z(i)-sig(i)) + sig(i);%���漰��ƽ�Ƶ�ÿһ����+1,-1
end
rate = sum(H(P));
dist = 0.5*double(sum(H.*(C.^2+(C+PK).^2))); 