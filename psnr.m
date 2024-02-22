function [PSNR, MSE] = psnr(X, Y)
% 计算峰值信噪比PSNR、均方根误差MSE
% 如果输入Y为空，则视为X与其本身来计算PSNR、MSE
%输入：
%   X，Y为两幅图像
%输出：
%   PSNR，峰值信噪比
%   MSE，均方误差
if nargin<2
    error('Not enough input arguments!');
else
    if any(size(X)~=size(Y))
        error('The input size is not equal to each other!');
    end
    D = int16(X)-int16(Y);
end
MSE = sum(D(:).*D(:))/numel(X);
PSNR = 10*log10(255^2/MSE);