function [PSNR, MSE] = psnr(X, Y)
% �����ֵ�����PSNR�����������MSE
% �������YΪ�գ�����ΪX���䱾��������PSNR��MSE
%���룺
%   X��YΪ����ͼ��
%�����
%   PSNR����ֵ�����
%   MSE���������
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