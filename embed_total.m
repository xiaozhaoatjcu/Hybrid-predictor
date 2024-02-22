  function [result,stego_img] = embed_total(img,pred_img,pl,opl,method,side,ds,dm,sect_switch,auxiliary,x_Edge)
%结合X/O集合像素均值预测，预测差值排序后分割，选点的综合实现
%输入：
%   img，原始图像
%   half_stego_img, 嵌入了一半秘密信息的图像
%   pl，总隐秘信息
%   opl, 剩余待嵌入信息
%   method，使用的排序方法，0为多特征，1为Sachnev特征，2为Hwang特征
%   twiceL二次排序
%输出：
%   result，最终结果和中间变量的记录
%   stego_img，载秘图象
%   PSNR，峰值信噪比


stego_img = [];
PSNR = 0;
                                                                                             %检测输入 
if numel(pl) > floor(numel(img)*0.9)                                       %限制最大嵌入率为0.9bpp
    error('Too much payload');
end


[r,c] = size(img);                                                         %检测图像尺寸，限制最小尺寸及偶数长宽
if r<16  || c<16 || mod(r,2) || mod(c,2)
    error('Invalid image size');
end

result.ok = 0;


%在O集合中嵌入一半隐秘信息
[o_result] = embed_in_one_set(img,pred_img, opl, 1, method,side,ds,dm,sect_switch);
if ~o_result.ok
    stego_img = zeros(r,c);
    stego_img = uint8(stego_img) ; 
    return;
end
%%%%%%
%%处理边缘信息
Edge_finite=numel(x_Edge)+auxiliary;
o_result.Edge=[o_result.PZ_infor o_result.side_info  o_result.LB];%%O系列的边缘信息
if numel(o_result.Edge)<=Edge_finite
%     fprintf('X边缘信息过大');
%     stego_img = zeros(r,c);
%     stego_img = uint8(stego_img) ;
% error('O系列边缘信息过大，错误位于embed主程序之中');
% end
    tmpdata=1;
    for tmp_r=1:r
         for tmp_c=1:c
             if (tmp_r>=side+1)&&(tmp_r<=r-side)&&(tmp_c>=side+1)&&(tmp_c<=c-side)
                 continue
             end
             if mod(tmp_r+tmp_c,2)==1
             %%奇数为O系列
            o_result.stego_img(tmp_r,tmp_c)=bitset(o_result.stego_img(tmp_r,tmp_c),1,o_result.Edge(tmpdata));
             tmpdata=tmpdata+1;
             else
             end
             if tmpdata>numel(o_result.Edge)
             break
             end
         end
         if tmpdata>numel(o_result.Edge)
         break
         end
    end
    %%处理边缘信息
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    o_result.dist = sum(sum((double(o_result.stego_img) - ...                  %O集合失真
        double(img)).^2));
    stego_img = uint8(o_result.stego_img);

    [PSNR, ~] = psnr(img, stego_img);                                        %计算峰值信噪比并记录中间变量

    result.O = o_result;
    result.ok = 1;
    result.PSNR = PSNR;
else
    stego_img = zeros(r,c);
    stego_img = uint8(stego_img) ; 

end
%%
function [result] = embed_in_one_set(img,pred_img, pl, set, method,side,ds,dm,sect_switch)
%在X/O集合中嵌入
%输入：
%   img，图像
%   pl，隐秘信息
%   set，集合标志，0为X集合，1为O集合
%   method，使用的排序方法，0为10个特征，1为Sachnev特征，2为Hwang特征，3为Li特征，4为30个特征
%输出
%   side特征选取的范围（3*3，5*5，7*7）
%   result，嵌入后的结果和中间变量

ok = 0;                                                                    %能否嵌入的标志
rate = numel(pl);                                                       

[pred_errors_cnn, pred_img_cnn, err_pos_cnn] = predict_one_set(img, pred_img,set,side);              %计算预测值和预测差值
Hi_cnn=hist(pred_errors_cnn,-255:255);
% His_cnn=histogram(pred_errors_cnn,-255:255);% 画出cnn预测误差直方图
tt1=clock();
[sorted_errors, sorted_err_pos,sort_side_info,record_info,cmT,sortsurS,corr,featurenumber,pred_img,pred_errors,err_pos] = ...       %预测差值排序
    sort_one_set(img,pred_img_cnn,pred_errors_cnn, err_pos_cnn ,set,method,side,dm);
Hi=hist(pred_errors,-255:255);
% His=histogram(pred_errors,-255:255);% 画出最终预测误差直方图
tt2=clock();
fea_T = etime(tt2,tt1);
tt1 = clock();
if ds==1
    time1=clock();
    if sect_switch==0
     [H,L,Length]=per_cut_histogram(sorted_errors,sortsurS);
    else
     [H,L,Length]=same_cut_histogram(sorted_errors);   
    end
    [Map_PZ,success,Length,MP,BIGDATA,S_Rate,S_dist]=LiSelect(rate,H,L,Length);%%Li的穷举法
    if success == 0      %如果不能搜索到满足要求的点，则返回
        result.ok = ok;
    return;
    end
    time2=clock();
    pz_time=etime(time2,time1);
    result.S_Rate=S_Rate;
    result.S_dist=S_dist;
    result.pz_time=pz_time;
    result.BIGDATA=BIGDATA;                                              %穷举法所用的情况
    %%%%
    %%增加边缘信息
    a=size(Map_PZ,1);
    PZ_infor= toBitSeq(a,5);%%Li的选点使用了16个直方图故只要存放五位边缘信息
    PZ_layer=sum(Map_PZ(:,end));%%总共使用了多少对点
    PZ_layer=toBitSeq(PZ_layer,10);%%用了多少对峰值点和零值点
    PZ_infor=[PZ_infor PZ_layer];%%前面五位多少个直方图（默认为16个）6-15（用了多少对点）
    b=sum(Map_PZ(:,end)==0);%%多少个直方图没有被使用
    a=a-b;%%有多少个直方图被使用，Li的方法直方图都是连续被使用的
    for i=1:a
    tmp_P=Map_PZ(i,1);
        if tmp_P<0
        PZ_infor=[PZ_infor 1];
        else
        PZ_infor=[PZ_infor 0];    
        end
    tmp_P=toBitSeq(abs(tmp_P),8);
    PZ_infor=[PZ_infor tmp_P];
    tmp_P=Map_PZ(i,2);
        if tmp_P<0
        PZ_infor=[PZ_infor 1];
        else
        PZ_infor=[PZ_infor 0];    
        end
    tmp_P=toBitSeq(abs(tmp_P),8);
    PZ_infor=[PZ_infor tmp_P];
    end
    result.PZ_infor=PZ_infor;
elseif ds==3
    [length,layer,P,Z,idx,Dist,DS,minPer,calCount] = devideAndPick2(sorted_errors,rate,99);
    success=1;
        if length==1
            success=0;
        end
        if success == 0      %如果不能搜索到满足要求的点，则返回
            result.ok = ok;
        return;
        end
    Length=[0 length numel(sorted_errors)];
    Map_PZ=zeros(2,2*numel(P)+1);
    Map_PZ(1,1:numel(P))=P;
    Map_PZ(1,numel(P)+1:2*numel(P))=Z;
    Map_PZ(1,2*numel(P)+1)=numel(P);
    Map_PZ(2,1:numel(P))=255;
    Map_PZ(2,numel(P)+1:2*numel(P))=255;
    Map_PZ(2,2*numel(P)+1)=0;
    [MP]= getModifyMap(P,Z);
    tmpMP=zeros(511,2);
        for i = -255:255
        tmpMP(i+256,:) = tmpMP(i+256,:) + i;                                         
        end
    MP=[MP tmpMP];
    result.layer=layer;%%峰值点和零值点的个数
    result.per=idx;%%最优百分比
    result.Dist=Dist;%%最优分割点处失真
    result.DS=DS;%%全部分割点处失真
    result.minPer=minPer;%%最小失真百分比
    result.calCount=calCount;%%测试百分比的个数
%%%%
%%存放边缘信息
PZ_infor=toBitSeq(length,20);
        if length>2^20
        errors('截点位置信息过长');
        end
PZ_layer=sum(Map_PZ(:,end));%%总共使用了多少对点
PZ_layer=toBitSeq(PZ_layer,10);%%用了多少对峰值点和零值点
PZ_infor=[PZ_infor PZ_layer];%%前面五位多少个直方图（默认为16个）6-15（用了多少对点）
tmpdata=Map_PZ(1,end);%%记录所有的峰值点和零值点
    for i=1:2*tmpdata
        tmpdata1=Map_PZ(1,i);
        if tmpdata1<0
        PZ_infor=[PZ_infor 1];
        else
        PZ_infor=[PZ_infor 0];
        end
        tmpdata1=toBitSeq(abs(tmpdata1),8);
        PZ_infor=[PZ_infor tmpdata1];
    end
result.PZ_infor=PZ_infor;
%%存放边缘信息
%%%%

else
    error('ds选点错误');
end
tt2 = clock();
[stego_img,LB,ok] = real_embed(sorted_errors, pred_img, img, Length,sorted_err_pos,...
   pl,Map_PZ,MP);
if ok==0
    error('嵌入信息real_embed秘密信息没有嵌入完成');
end

%%计算使用峰值点和零值点的总数
n=size(Map_PZ,2);
sumlayer=sum(Map_PZ(:,n));


perT = etime(tt2,tt1);
[H_definite_infor]=simulate_entropy(Length,sorted_errors);
result.H_definite_infor=H_definite_infor;
%记录结果和中间变量
result.stego_img = stego_img;                                              %载秘图像
result.ok = ok;                                                            %嵌入成功标志
result.Map_PZ=Map_PZ;                                            %%存放峰值点和零值点
result.MP=MP;                                                           %%累计平移量
result.sumlayer=sumlayer;                                          %%总共使用的点对数
result.side_info = sort_side_info;                                         %排序用的边信息（特征权值）
result.record_info = record_info;                                          %排序时的中间变量
result.rate=rate;                                                          %嵌入容量
result.perT = perT;                                                        %搜索分割百分比的时间
result.cmT = cmT;                                                          %计算邻域复杂度测量值的时间
result.fea_T=fea_T;                                                                %排序的时间
result.sorted_errors = sorted_errors;                                      %排序后的预测差值
result.sortsurS=sortsurS;                                              %记录排序预测误差与其相对应排序后的像素的各种特征值
result.corr=corr;                                                           %未排序预测误差的绝对值与多特征值的相关系数
result.Length=Length;                                                    %截取拟合特征的长度
result.LB=LB;                                                                 %生成的位图
result.featurenumber=featurenumber;                           %自动生成特征的个数
%%
function [stego_img,LB,ok] = real_embed(sorted_errors, pred_img, cover_img, Length, sorted_err_pos,...
     pl,Map_PZ,MP)
%%输入
%%sorted_errors排序的预测误差
%%pred_img预测图像
%%cover_img原图像
%%Length截段的长度
%%sorted_err_pos排序预测误差的坐标
%%pl嵌入信息
%%Map_PZ截断的PZ点
%%MP每个截断点的累计平移量

stego_img = double(cover_img);
pred_img = double(pred_img);                             
m=numel(pl);
n=numel(Length)-1;
d=size(Map_PZ,2);
j=0;
LB=zeros(1,n);
for i=1:n
 tmp_LB=[];
 tmp_sorted_errors=sorted_errors(Length(i)+1:Length(i+1));   
 tmp_sorted_err_pos=sorted_err_pos(Length(i)+1:Length(i+1),:);           
 mp=MP(:,2*i-1:2*i);
 tmplayer=Map_PZ(i,d);
 if tmplayer==0
     continue
 end
 P=Map_PZ(i,1:tmplayer);
 d1=numel(tmp_sorted_errors);
      for s=1:d1 
              key=ismember(tmp_sorted_errors(s),P);
                 if j<m%%嵌入容量是否超过pl容量
                    if  key%%是否为峰值点
                        j=j+1;
                       fake=pred_img(tmp_sorted_err_pos(s,1),tmp_sorted_err_pos(s,2))+mp(tmp_sorted_errors(s)+256,1+pl(j));
                    else
                       fake=pred_img(tmp_sorted_err_pos(s,1),tmp_sorted_err_pos(s,2))+mp(tmp_sorted_errors(s)+256,1);
                    end
                 else
                       fake=pred_img(tmp_sorted_err_pos(s,1),tmp_sorted_err_pos(s,2))+mp(tmp_sorted_errors(s)+256,1); 
                 end
                
                if fake<=tmplayer-1%%防止溢出同时生成位图LB
                    if fake<0;
                        tmp_LB=[tmp_LB 1];
                       fake=fake+tmplayer;
                    else
                        tmp_LB=[tmp_LB 0];
                    end
                end
            
                if fake>=255-tmplayer+1
                     if fake>255;
                         tmp_LB=[tmp_LB 1];
                         fake=fake-tmplayer;
                     else
                         tmp_LB=[tmp_LB 0];
                     end
                end  
                stego_img(tmp_sorted_err_pos(s,1),tmp_sorted_err_pos(s,2))=fake;
      end
        if sum(tmp_LB==1)>0
            LB(i)=1;
            [LB]=real_embed_LB_compose(LB,tmp_LB);
        else
            LB(i)=0;
        end
end

if j<m
    ok=0;
else
    ok=1;
end
%%
function[LB]=real_embed_LB_compose(LB,tmp_LB)
tmp_data=sum(tmp_LB==1);
m=numel(tmp_LB);
m=dec2bin(m);
m=numel(m);
if tmp_data>8
LB=[LB 0 tmp_LB];%%不进行压缩
else
[~,pos]=find(tmp_LB==1);
n=numel(pos);
n=toBitSeq(n-1,3);
LB=[LB 1 n];
    for i=1:numel(pos);
    T= toBitSeq(pos(i),m);
    LB=[LB T];
    end
end
%%
function [si] = get_embed_side_info(P,Z,length)
%编码边信息
%输入：
%   P，峰值点
%   Z，零点
%   length，实际使用差值数
%输出：
%   si，编码后边信息，包括峰值点数，峰值点-零点，实际使用差值数

si = zeros(1,numel(P)*18 + 17);
for i = 1:numel(P)
    si(1,i*18-16:i*18-9) = toBitSeq(P(i),8);
    si(1,i*18-7:i*18) = toBitSeq(Z(i),8);
    if P(i)<0
        si(1,i*18-17) = 1;
    end
    if Z(i)<0
        si(1,i*18-8) = 1;
    end
end
si(1,numel(P)*18+1:numel(si)) = toBitSeq(length,17);
si = [toBitSeq(numel(P),5),si];
%%
function [mp] = getModifyMap(P,Z)
%根据每一层平移和嵌入的叠加构造差值映射表
%输入：
%   P，峰值点
%   Z，零点
%输出：
%   mp，差值映射表，若差值i不是峰值点，则mp（i+256，1） = mp（i+256，2）为平移后差值,
%       否则mp（i+256，1）为嵌入0后的差值，mp（i+256，2）为嵌入1后的差值

mp = zeros(511,2);
P = P+256;
Z = Z+256;
for i=1:numel(P)
    if P(i) < Z(i)
        if P(i)+1 <= Z(i)-1
            mp(P(i)+1:Z(i)-1,:) = mp(P(i)+1:Z(i)-1,:) + 1;                 %平移
        end
        mp(P(i),2) = mp(P(i),2) + 1;                                       %嵌入
    else
        if P(i)-1 >= Z(i)+1
            mp(Z(i)+1:P(i)-1,:) = mp(Z(i)+1:P(i)-1,:) - 1;                 %平移
        end
        mp(P(i),2) = mp(P(i),2) - 1;                                       %嵌入
    end
end
for i = -255:255
    mp(i+256,:) = mp(i+256,:) + i;                                         %预测误差加上修改量后的预测误差
end


%%
function [pred_errors_cnn, pred_img_cnn, err_pos_cnn] = predict_one_set(img, pred_img,set,side)
%获得预测图像和预测差值
%输入：
%   img：原图像
%   set：X/O集，0代表X集
%   side特征选取的范围（3*3，5*5，7*7）
%输出：
%   pred_errors：预测差值数组
%   pred_img：预测图像
%   img_pos：预测差值对应的图像位置

[r,c] = size(img);
img = double(img);
pre_img = double(pred_img);
pred_img_cnn = img;
pred_img_cnn(side+1:r-side,side+1:c-side) = pre_img(side+1:r-side,side+1:c-side);
% 将原图像中另一个集合的像素值与预测的集合组合成一个集合的预测图
for row = side+1:r-side
     for col = side+1:c-side
        if mod(row + col, 2) == 0
            pred_img_cnn(row,col) = img(row,col);
         end
     end
end


pred_errors_cnn = zeros(1,(r-side*2)*(c-side*2)/2);
err_pos_cnn = zeros(numel(pred_errors_cnn),2);

k = 1;
for row = side+1:r-side
    for col = side+1:c-side
        if mod(row + col, 2) == set
            pred_errors_cnn(k) = img(row,col) - pred_img_cnn(row,col);
            err_pos_cnn(k,1) = row;
            err_pos_cnn(k,2) = col;
            k = k + 1;
        end
    end
end

%%
function [sorted_errors, sorted_err_pos,sort_side_info,record_info,cmt,sortsurS,corr,featurenumber,pred_img,pred_errors,err_pos] = sort_one_set(img,pred_img_cnn,pred_errors, err_pos ,set,method,side,dm)
%选择最优特征权重，排序预测差值
%输入：
%   img，预测值图像
%   pred_errors，预测差值序列
%   err_pos，预测差值在图像中位置
%   set，X/O集标志，0代表X集
%   method,采用的排序特征，0代表多特征，1代表Sachnev的特征，2代表Hwang的特征
%输出：
%   sorted_errors，排序后的预测差值序列
%   sorted_err_pos，排序后差值在图像中位置
%   sort_side_info，二进制序列，排序的优化参数
%   record_info，记录用于观察的信息
%   cmt,计算特征值所用的时间
%   sortsurS,排序后的预测误差与排序后的像素与其相对应的特征值
%   corr，未排序的预测误差的绝对值与多特征值的相关系数
%   featurenumber特征个数

%初始化
[r,c] = size(img);
img = double(img);
sorted_errors = pred_errors;
sorted_err_pos = err_pos;


if method==0
        [surS,w,cmt,sur,ms,originw]=old_feature(img,pred_errors,r,c,set,side,dm);
        %%tmp_w权重值
        %%tmp_sur未排序的预测误差
        %%tmp_originw未整数化的权重值
        %%ms值的数量特征
        %%tmp_cmt所用的时间
        %%surS所有的单特征值
        sortsurS=zeros(ms+2,numel(sur));
         [~,pos] = sort(sur);
             for j = 1:numel(sur)
                 sorted_errors(1,j) = pred_errors(1,pos(j));
                 sorted_err_pos(j,1) = err_pos(pos(j),1);
                 sorted_err_pos(j,2) = err_pos(pos(j),2);
                 sortsurS(1:ms,j)=surS(1:ms,pos(j));   
                 sortsurS(ms+1,j)=sur(1,pos(j));%%排序后的多特征值
             end   
        featurenumber=10;
        corr=[];
        record_info.w = w;
        record_info.originw=originw;
        sort_side_info = getWSeq(w);
        sortsurS(ms+2,:)=sorted_errors;
        absErr = abs(pred_errors);
        record_info.corr = corrcoef(absErr,double(w)*surS);
        corr=record_info.corr;
elseif ismember(method,[1 2 3])
    [surS,cmt,sur]=single_feature(img,pred_errors,r,c,set,method,side);
    originw=[];
    featurenumber=1;
    w=1;
    ms=1;
else
    error('method选点错误');
end

if ismember(method,[1 2 3])
%排序
[~,pos] = sort(sur);
sortsurS=zeros(ms+2,numel(sur));

% 根据复杂度大小选择不同预测器，菱形预测用0表示，cnn预测器用1表示
[pred_errors, pred_img, err_pos, prediction_map] = prediction_select(img,pred_img_cnn,side,set,method);
record_info.prediction_map = prediction_map;

    for i = 1:numel(sur)
      sorted_errors(1,i) = pred_errors(1,pos(i));
      sorted_err_pos(i,1) = err_pos(pos(i),1);
      sorted_err_pos(i,2) = err_pos(pos(i),2);
      sortsurS(1:ms,i)=surS(1:ms,pos(i));   
      sortsurS(ms+1,i)=sur(1,pos(i));%%排序后的多特征值
    end
sortsurS(ms+2,:)=sorted_errors;

sort_side_info = getWSeq(w);
record_info.w = w;
record_info.originw=originw;

absErr = abs(pred_errors);
record_info.corr = corrcoef(absErr,double(w)*surS);
corr=record_info.corr;
end

%%
function[surS,cmt,sur]=single_feature(img,pred_errors,r,c,set,method,side)
%%输入
%%img输入图像
%%pred_errors预测误差
%%r,c图像的长和宽
%%set，0表示x集合，1表示o集合
%%method 1表示Sachnev3*3差值绝对值的方差
%%method2 表示Hwang3*3方差
%%method3表示Li的特征
%%side表示选用特征的宽度，指图像边上side长度的像素不执行嵌入
t1=clock();
sur=zeros(1,numel(pred_errors));
sur1=sur;
k=1;
for row = side+1:r-side
    for col = side+1:c-side
        if mod(row+col,2) == set
            if method==1
            surround = [img(row,col-1),img(row-1,col),img(row,col+1),img(row+1,col)];%3*3邻域像素
            dist_near = abs(diff([surround,img(row,col-1)]));%邻域像素的差值
            sur1(1,k) = var(dist_near,1);
            end
            if method==2
            surround = [img(row,col-1),img(row-1,col),img(row,col+1),img(row+1,col)];%3*3邻域像素
            sur1(1,k)= var(surround,1);
            end
            if method==3
                surround3=[abs(img(row,col-1)-img(row+1,col-1)),abs(img(row+1,col-1)-img(row+2,col-1)),abs(img(row+1,col)-img(row+2,col))...
             ,abs(img(row,col+1)-img(row+1,col+1)),abs(img(row+1,col+1)-img(row+2,col+1)),abs(img(row-1,col+2)-img(row,col+2))...
             ,abs(img(row,col+2)-img(row+1,col+2)),abs(img(row+1,col+2)-img(row+2,col+2)),abs(img(row,col+1)-img(row,col+2))...
             ,abs(img(row+1,col-1)-img(row+1,col)),abs(img(row+1,col)-img(row+1,col+1)),abs(img(row+1,col+1)-img(row+1,col+2))...
             ,abs(img(row+2,col-1)-img(row+2,col)),abs(img(row+2,col)-img(row+2,col+1)),abs(img(row+2,col+1)-img(row+2,col+2))];
            sur1(1,k)=sum(surround3);
            end
            k=k+1;
        end
    end
end
% k=1;
% for row = side+1:r-side                           %%使生成的特征值不再重复
%     for col = side+1:c-side
%         if mod(row+col,2)==set
%     sur(k)=sur(k)+(row+(col-1)*r)/10^(15);
%     k=k+1;
%         end
%     end
% end
sur=sur1;
surS=sur1;
t2=clock();
cmt=etime(t2,t1);

%%
function[surS,w,cmt,sur,ms,originw]=old_feature(img,pred_errors,r,c,set,side,dm)
%特征值序列
sur = zeros(1,numel(pred_errors));
%12个特征值序列
sur1 = sur;
sur2 = sur;
sur5 = sur;
sur6 = sur;
sur7 = sur;
sur8 = sur;
sur9 = sur;
sur10 = sur;
sur11 = sur;
sur12 = sur;

%对角的二阶差值算子
filterBlk = [1,0,1;0,-4,0;1,0,1];
%用于计算梯度的像素的位置相对中间像素的偏差
gradOffset = [-1,0,1,0;-2,1,0,1;-1,2,1,2;0,1,2,1];
gradOffset2 = [-1,-2,1,-2;-2,-1,0,-1;-1,0,1,0;0,-1,2,-1];
%梯度差值
gradDiff = zeros(6,2);
gradDiff2 = zeros(6,2);
k = 1;
tt1 = clock();
for row = side+1:r-side
    for col = side+1:c-side
        if mod(row+col,2) == set
            %3*3邻域像素
            surround = [img(row,col-1),img(row-1,col),img(row,col+1),img(row+1,col)];
            %邻域像素的差值
            dist_near = abs(diff([surround,img(row,col-1)]));
                %计算各个相邻像素的梯度
                grad = [img(row+gradOffset(1,1),col+gradOffset(1,2))-surround(1),img(row+gradOffset(1,3),col+gradOffset(1,4))-surround(1);
                    img(row+gradOffset(2,1),col+gradOffset(2,2))-surround(2),img(row+gradOffset(2,3),col+gradOffset(2,4))-surround(2);
                    img(row+gradOffset(3,1),col+gradOffset(3,2))-surround(3),img(row+gradOffset(3,3),col+gradOffset(3,4))-surround(3);
                    img(row+gradOffset(4,1),col+gradOffset(4,2))-surround(4),img(row+gradOffset(4,3),col+gradOffset(4,4))-surround(4)];
                grad2 = [img(row+gradOffset2(1,1),col+gradOffset2(1,2))-surround(1),img(row+gradOffset2(1,3),col+gradOffset2(1,4))-surround(1);
                    img(row+gradOffset2(2,1),col+gradOffset2(2,2))-surround(2),img(row+gradOffset2(2,3),col+gradOffset2(2,4))-surround(2);
                    img(row+gradOffset2(3,1),col+gradOffset2(3,2))-surround(3),img(row+gradOffset2(3,3),col+gradOffset2(3,4))-surround(3);
                    img(row+gradOffset2(4,1),col+gradOffset2(4,2))-surround(4),img(row+gradOffset2(4,3),col+gradOffset2(4,4))-surround(4)];
                %计算各个相邻像素的梯度差
                gk = 1;
                for p = 1:3
                    for q = p+1:4
                        gradDiff(gk,:) = grad(p,:) - grad(q,:);
                        gradDiff2(gk,:) = grad2(p,:) - grad2(q,:);
                        gk = gk+1;
                    end
                end
            

            %Hwang和Sachnev所用的特征值
            sur6(1,k) = var(surround,1);%Hwang
            sur7(1,k) = var(dist_near,1);%Sachnev
            sur1(1,k) = max(surround);
            sur2(1,k) = min(surround);
            sur5(1,k) = sum(dist_near);
            sur8(1,k) = abs(sum(sum(img(row-2:row,col-1:col+1).*filterBlk))) + abs(sum(sum(img(row:row+2,col-1:col+1).*filterBlk))) + abs(sum(sum(img(row-1:row+1,col-2:col).*filterBlk))) + abs(sum(sum(img(row-1:row+1,col:col+2).*filterBlk)));
            sur9(1,k) = sum(sum(gradDiff.^2,2).^0.5);
            sur10(1,k) = sum(sum(gradDiff2.^2,2).^0.5);
            sur11(1,k) = sum(sum(grad.^2,2).^0.5);
            sur12(1,k) = sum(sum(grad2.^2,2).^0.5);
            
            k = k + 1;
        end
    end
end

%将所有特征值缩放到[0,1000]的范围内
[sur6,~] = mapminmax(sur6);
sur6 = (sur6+1)*500;
[sur7,~] = mapminmax(sur7);
sur7 = (sur7+1)*500;
[sur1,~] = mapminmax(sur1);
sur1 = (sur1+1)*500;
[sur2,~] = mapminmax(sur2);
sur2 = (sur2+1)*500;
[sur5,~] = mapminmax(sur5);
sur5 = (sur5+1)*500;
[sur8,~] = mapminmax(sur8);
sur8 = (sur8+1)*500;
[sur9,~] = mapminmax(sur9);
sur9 = (sur9+1)*500;
[sur10,~] = mapminmax(sur10);
 sur10 = (sur10+1)*500;
 [sur11,~] = mapminmax(sur11);
 sur11 = (sur11+1)*500;
 [sur12,~] = mapminmax(sur12);
 sur12 = (sur12+1)*500;

surS = [sur1;sur2;sur5;sur6;sur7;sur8;sur9;sur10;sur11;sur12];
if dm
    absErr = abs(pred_errors);%%计算权重
    errNum = numel(absErr);
    lb=ones(1,10).*(-255);%%权重上下限
    ub=ones(1,10).*(255);
    options=optimset('Algorithm','active-set');
    w = fmincon(@(w) calCorr(w,absErr,surS,errNum),ones(1,10),[],[],[],[],lb,ub,[],options);
else
    V=surS';
    I=eye(size(V,2));
    y=(V'*V+0.00001*I)^(-1)*V'*abs(pred_errors');
    w=y';
end
originw=w;
w = int32(w/max(abs(w))*255);%范围为-255~255
% if di
% sur=  double(originw)*surS;
% intergersur=double(w)*surS;
% else
sur = double(w)*surS;
% intergersur=double(w)*surS;
% end


k=1;
for row = side+1:r-side                           %%使生成的特征值不再重复
    for col = side+1:c-side
        if mod(row+col,2)==set
    sur(k)=sur(k)+(row+(col-1)*r)/10^(15);
    k=k+1;
        end
    end
end
tt2 = clock();
cmt = etime(tt2,tt1);
ms=10;%%our方法所使用的特征
%%
function [res] = calCorr(w,absErr,surS,k)
%计算权重为w，特征值为surS时和预测差值的相关系数
%输入：
%   w，权重
%   err，差值序列
%   surS，全部特征值
%   k，差值个数
%输出：
%   res，相关系数
fea = double(w)*surS;
meanF = mean(fea);
res  = 0-(absErr*fea'/k - mean(absErr)*meanF)/(sqrt(mean((fea).^2) - (meanF)^2)*sqrt(mean((absErr).^2)-(mean(absErr)).^2));
%%
function [BS]  = toBitSeq(x,l)%将输入的x即w(i)转换为一个l位(即8位)二进制数
%由x得到x绝对值的二进制序列
%输入:
%   x，整数
%   l，序列长度
%输出:
%   BS，二进制序列
BS = zeros(1,l);
x = uint32(abs(x));
for i=1:l
    BS(1,i) = bitget(x,l-i+1);%将x转换为二进制数，取低位开始的i位数
end
%%
function [BS] =  getWSeq(w)
%将权重编码为01序列
%输入：
%   w，权重
%   BS，编码后序列
numW = numel(w);
BS = zeros(1,numW*9);
for i = 1:numW
    BS(1,i*9-7:i*9) = toBitSeq(w(i),8);
    if w(i) < 0
        BS(1,i*9-8) = 1;
    else
        BS(1,i*9-8) = 0;
    end
end




