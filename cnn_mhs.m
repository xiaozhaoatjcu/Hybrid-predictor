function stego_img = cnn_mhs(img,pred_img,m,set)

%-----------����ֱ��ͼǶ��-----------

%���룺
%   img��ԭʼͼ��
%   pred_img������cnnģ��Ԥ��ͼ��
%   m��������Ϣ��m bpp��
%   set��ͼ���cross���ϻ���dot����
%�����
%   stego_img������Ƕ��ͼ��

time1=clock();
img = double(img);   %ԭʼͼ��
pred_img = double(pred_img);   %cnnԤ��ͼ��

auxiliary=250;%������Ϣ
pl = getPayload(m);   %��Ϣ     ��������Ϊ�ɵ��ء���ť��
method =3;    %���򷽷�0Ϊ10������1ΪSachnev������2ΪHwang������GA����3ΪLi����
sect_switch=0;%%1Ϊ�ȷ֣�0Ϊ�ȱ�����
side=2;%ѡȡ�����ķ�Χ 
ds=3;%ѡ��ķ�ʽ��1ΪLi�ķ��䷽ʽ����ֱ��ͼ����3Ϊѡȡ���Űٷֱȣ���ֱ��ͼ��
dm=1;%%(method=0) dm=1Ϊѡ��fmincon,dm=0Ϊѡ����Ƶ��ʽ 
 
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


