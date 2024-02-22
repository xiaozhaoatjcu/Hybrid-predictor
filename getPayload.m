function [pl] = getPayload(bpp,size)
%����Ƕ���ʺ�ͼ���С����ȡһ��01���ȷֲ���������Ϣ
%���룺
%   bpp��Ƕ���ʣ�����С��1
%   size��ͼ��߳�
%�����
%   pl��һ��01����
if nargin<2
    size=512;
end
if bpp>=1
    error('bpp is larger than 1');
end
pl = rand(1,uint32(size*size*bpp))>0.5;

m = sum(pl==0);  
n = sum(pl==1);
if n>m
    pl = ~pl;
end
