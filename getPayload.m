function [pl] = getPayload(bpp,size)
%根据嵌入率和图像大小，获取一个01均匀分布的隐秘信息
%输入：
%   bpp：嵌入率，必须小于1
%   size，图像边长
%输出：
%   pl，一个01数组
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
