function stego_img = cnn_mhs(img,pred_img,m,set)

%-----------进行直方图嵌入-----------

%输入：
%   img：原始图像
%   pred_img：利用cnn模型预测图像
%   m：秘密信息（m bpp）
%   set：图像的cross集合或者dot集合
%输出：
%   stego_img：最终嵌密图像

time1=clock();
img = double(img);   %原始图像
pred_img = double(pred_img);   %cnn预测图像

auxiliary=250;%辅助信息
pl = getPayload(m);   %信息     以下数据为可调控【按钮】
method =3;    %排序方法0为10特征，1为Sachnev特征，2为Hwang特征（GA），3为Li特征
sect_switch=0;%%1为等分；0为等比例分
side=2;%选取特征的范围 
ds=3;%选点的方式，1为Li的分配方式（多直方图），3为选取最优百分比（单直方图）
dm=1;%%(method=0) dm=1为选择fmincon,dm=0为选择视频公式 
 
if set==0
    [result,opl,stego_img,PSNR,x_Edge] = embed_half(img,pred_img,pl,method,side,ds,dm,sect_switch,auxiliary); 
    save('opl.mat', 'opl')
    save('x_Edge.mat', 'x_Edge')
    save('pl.mat', 'pl')
else
    load('D:\pengyi\master_student\work_CNN\xwcnnmhs_final\opl');
    load('D:\pengyi\master_student\work_CNN\xwcnnmhs_final\x_Edge');
    load('D:\pengyi\master_student\work_CNN\xwcnnmhs_final\pl');
    [result,stego_img] = embed_total(img,pred_img,pl,opl,method,side,ds,dm,sect_switch,auxiliary,x_Edge);
    imwrite(stego_img, 'D:\pengyi\master_student\work_CNN\xwcnnmhs_final\stego_images_dir\stego_img.bmp');
end
time2=clock();
time=etime(time2,time1);


