function [H] = diamond(f, xo)
%菱形预测
%输入：
%   f：图像，默认长宽都是偶数
%   xo：X/O集标志
%输出：
%   H：预测差值直方图

[r,c] = size(f);
L = int16(zeros(1,r*c/2));
tf = double(f);
k = 1;
if xo==1
    b = 0;
else
    b = 1;
end
for i=2:r-1
    for j=2:c-1
        if mod(i+j,2)==b
            sum = tf(i-1,j) + tf(i+1,j) + tf(i,j-1) + tf(i,j+1);
            L(k) = tf(i,j) - round(sum/4);
            k = k + 1;
        end
    end
    if mod(i+1,2)==b
        sum = tf(i-1,1) + tf(i+1,1) + tf(i,2);
        L(k) = tf(i,1) - round(sum/3);
        k = k + 1;
    end
    if mod(i+c,2)==b
        sum = tf(i-1,c) + tf(i+1,c) + tf(i,c-1);
        L(k) = tf(i,c) - round(sum/3);
        k = k + 1;
    end
end
for j=2:c-1
    if mod(j+1,2)==b
        sum=tf(1,j-1) + tf(1,j+1) + tf(2,j);
        L(k) = tf(1,j)-round(sum/3);
        k = k+1;
    end
    if mod(j+r,2)==b
        sum = tf(r,j-1) + tf(r,j+1) + tf(r-1,j);
        L(k) = tf(r,j) - round(sum/3);
        k = k + 1;
    end
end
if b==0
    L(k) = tf(1,1) - round((tf(1,2)+tf(2,1))/2);
    k = k + 1;
    L(k) = tf(r,c) - round((tf(r,c-1)+tf(r-1,c))/2);
else
    L(k) = tf(1,c) - round((tf(1,c-1)+tf(2,c))/2);
    k = k + 1;
    L(k) = tf(r,1) - round((tf(r,2)+tf(r-1,1))/2);
end
H=hist(L,-255:1:255);