  function [result,stego_img] = embed_total(img,pred_img,pl,opl,method,side,ds,dm,sect_switch,auxiliary,x_Edge)
%���X/O�������ؾ�ֵԤ�⣬Ԥ���ֵ�����ָѡ����ۺ�ʵ��
%���룺
%   img��ԭʼͼ��
%   half_stego_img, Ƕ����һ��������Ϣ��ͼ��
%   pl����������Ϣ
%   opl, ʣ���Ƕ����Ϣ
%   method��ʹ�õ����򷽷���0Ϊ��������1ΪSachnev������2ΪHwang����
%   twiceL��������
%�����
%   result�����ս�����м�����ļ�¼
%   stego_img������ͼ��
%   PSNR����ֵ�����


stego_img = [];
PSNR = 0;
                                                                                             %������� 
if numel(pl) > floor(numel(img)*0.9)                                       %�������Ƕ����Ϊ0.9bpp
    error('Too much payload');
end


[r,c] = size(img);                                                         %���ͼ��ߴ磬������С�ߴ缰ż������
if r<16  || c<16 || mod(r,2) || mod(c,2)
    error('Invalid image size');
end

result.ok = 0;


%��O������Ƕ��һ��������Ϣ
[o_result] = embed_in_one_set(img,pred_img, opl, 1, method,side,ds,dm,sect_switch);
if ~o_result.ok
    stego_img = zeros(r,c);
    stego_img = uint8(stego_img) ; 
    return;
end
%%%%%%
%%�����Ե��Ϣ
Edge_finite=numel(x_Edge)+auxiliary;
o_result.Edge=[o_result.PZ_infor o_result.side_info  o_result.LB];%%Oϵ�еı�Ե��Ϣ
if numel(o_result.Edge)<=Edge_finite
%     fprintf('X��Ե��Ϣ����');
%     stego_img = zeros(r,c);
%     stego_img = uint8(stego_img) ;
% error('Oϵ�б�Ե��Ϣ���󣬴���λ��embed������֮��');
% end
    tmpdata=1;
    for tmp_r=1:r
         for tmp_c=1:c
             if (tmp_r>=side+1)&&(tmp_r<=r-side)&&(tmp_c>=side+1)&&(tmp_c<=c-side)
                 continue
             end
             if mod(tmp_r+tmp_c,2)==1
             %%����ΪOϵ��
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
    %%�����Ե��Ϣ
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    o_result.dist = sum(sum((double(o_result.stego_img) - ...                  %O����ʧ��
        double(img)).^2));
    stego_img = uint8(o_result.stego_img);

    [PSNR, ~] = psnr(img, stego_img);                                        %�����ֵ����Ȳ���¼�м����

    result.O = o_result;
    result.ok = 1;
    result.PSNR = PSNR;
else
    stego_img = zeros(r,c);
    stego_img = uint8(stego_img) ; 

end
%%
function [result] = embed_in_one_set(img,pred_img, pl, set, method,side,ds,dm,sect_switch)
%��X/O������Ƕ��
%���룺
%   img��ͼ��
%   pl��������Ϣ
%   set�����ϱ�־��0ΪX���ϣ�1ΪO����
%   method��ʹ�õ����򷽷���0Ϊ10��������1ΪSachnev������2ΪHwang������3ΪLi������4Ϊ30������
%���
%   side����ѡȡ�ķ�Χ��3*3��5*5��7*7��
%   result��Ƕ���Ľ�����м����

ok = 0;                                                                    %�ܷ�Ƕ��ı�־
rate = numel(pl);                                                       

[pred_errors_cnn, pred_img_cnn, err_pos_cnn] = predict_one_set(img, pred_img,set,side);              %����Ԥ��ֵ��Ԥ���ֵ
Hi_cnn=hist(pred_errors_cnn,-255:255);
% His_cnn=histogram(pred_errors_cnn,-255:255);% ����cnnԤ�����ֱ��ͼ
tt1=clock();
[sorted_errors, sorted_err_pos,sort_side_info,record_info,cmT,sortsurS,corr,featurenumber,pred_img,pred_errors,err_pos] = ...       %Ԥ���ֵ����
    sort_one_set(img,pred_img_cnn,pred_errors_cnn, err_pos_cnn ,set,method,side,dm);
Hi=hist(pred_errors,-255:255);
% His=histogram(pred_errors,-255:255);% ��������Ԥ�����ֱ��ͼ
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
    [Map_PZ,success,Length,MP,BIGDATA,S_Rate,S_dist]=LiSelect(rate,H,L,Length);%%Li����ٷ�
    if success == 0      %�����������������Ҫ��ĵ㣬�򷵻�
        result.ok = ok;
    return;
    end
    time2=clock();
    pz_time=etime(time2,time1);
    result.S_Rate=S_Rate;
    result.S_dist=S_dist;
    result.pz_time=pz_time;
    result.BIGDATA=BIGDATA;                                              %��ٷ����õ����
    %%%%
    %%���ӱ�Ե��Ϣ
    a=size(Map_PZ,1);
    PZ_infor= toBitSeq(a,5);%%Li��ѡ��ʹ����16��ֱ��ͼ��ֻҪ�����λ��Ե��Ϣ
    PZ_layer=sum(Map_PZ(:,end));%%�ܹ�ʹ���˶��ٶԵ�
    PZ_layer=toBitSeq(PZ_layer,10);%%���˶��ٶԷ�ֵ�����ֵ��
    PZ_infor=[PZ_infor PZ_layer];%%ǰ����λ���ٸ�ֱ��ͼ��Ĭ��Ϊ16����6-15�����˶��ٶԵ㣩
    b=sum(Map_PZ(:,end)==0);%%���ٸ�ֱ��ͼû�б�ʹ��
    a=a-b;%%�ж��ٸ�ֱ��ͼ��ʹ�ã�Li�ķ���ֱ��ͼ����������ʹ�õ�
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
        if success == 0      %�����������������Ҫ��ĵ㣬�򷵻�
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
    result.layer=layer;%%��ֵ�����ֵ��ĸ���
    result.per=idx;%%���Űٷֱ�
    result.Dist=Dist;%%���ŷָ�㴦ʧ��
    result.DS=DS;%%ȫ���ָ�㴦ʧ��
    result.minPer=minPer;%%��Сʧ��ٷֱ�
    result.calCount=calCount;%%���԰ٷֱȵĸ���
%%%%
%%��ű�Ե��Ϣ
PZ_infor=toBitSeq(length,20);
        if length>2^20
        errors('�ص�λ����Ϣ����');
        end
PZ_layer=sum(Map_PZ(:,end));%%�ܹ�ʹ���˶��ٶԵ�
PZ_layer=toBitSeq(PZ_layer,10);%%���˶��ٶԷ�ֵ�����ֵ��
PZ_infor=[PZ_infor PZ_layer];%%ǰ����λ���ٸ�ֱ��ͼ��Ĭ��Ϊ16����6-15�����˶��ٶԵ㣩
tmpdata=Map_PZ(1,end);%%��¼���еķ�ֵ�����ֵ��
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
%%��ű�Ե��Ϣ
%%%%

else
    error('dsѡ�����');
end
tt2 = clock();
[stego_img,LB,ok] = real_embed(sorted_errors, pred_img, img, Length,sorted_err_pos,...
   pl,Map_PZ,MP);
if ok==0
    error('Ƕ����Ϣreal_embed������Ϣû��Ƕ�����');
end

%%����ʹ�÷�ֵ�����ֵ�������
n=size(Map_PZ,2);
sumlayer=sum(Map_PZ(:,n));


perT = etime(tt2,tt1);
[H_definite_infor]=simulate_entropy(Length,sorted_errors);
result.H_definite_infor=H_definite_infor;
%��¼������м����
result.stego_img = stego_img;                                              %����ͼ��
result.ok = ok;                                                            %Ƕ��ɹ���־
result.Map_PZ=Map_PZ;                                            %%��ŷ�ֵ�����ֵ��
result.MP=MP;                                                           %%�ۼ�ƽ����
result.sumlayer=sumlayer;                                          %%�ܹ�ʹ�õĵ����
result.side_info = sort_side_info;                                         %�����õı���Ϣ������Ȩֵ��
result.record_info = record_info;                                          %����ʱ���м����
result.rate=rate;                                                          %Ƕ������
result.perT = perT;                                                        %�����ָ�ٷֱȵ�ʱ��
result.cmT = cmT;                                                          %���������ӶȲ���ֵ��ʱ��
result.fea_T=fea_T;                                                                %�����ʱ��
result.sorted_errors = sorted_errors;                                      %������Ԥ���ֵ
result.sortsurS=sortsurS;                                              %��¼����Ԥ������������Ӧ���������صĸ�������ֵ
result.corr=corr;                                                           %δ����Ԥ�����ľ���ֵ�������ֵ�����ϵ��
result.Length=Length;                                                    %��ȡ��������ĳ���
result.LB=LB;                                                                 %���ɵ�λͼ
result.featurenumber=featurenumber;                           %�Զ����������ĸ���
%%
function [stego_img,LB,ok] = real_embed(sorted_errors, pred_img, cover_img, Length, sorted_err_pos,...
     pl,Map_PZ,MP)
%%����
%%sorted_errors�����Ԥ�����
%%pred_imgԤ��ͼ��
%%cover_imgԭͼ��
%%Length�ضεĳ���
%%sorted_err_pos����Ԥ����������
%%plǶ����Ϣ
%%Map_PZ�ضϵ�PZ��
%%MPÿ���ضϵ���ۼ�ƽ����

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
                 if j<m%%Ƕ�������Ƿ񳬹�pl����
                    if  key%%�Ƿ�Ϊ��ֵ��
                        j=j+1;
                       fake=pred_img(tmp_sorted_err_pos(s,1),tmp_sorted_err_pos(s,2))+mp(tmp_sorted_errors(s)+256,1+pl(j));
                    else
                       fake=pred_img(tmp_sorted_err_pos(s,1),tmp_sorted_err_pos(s,2))+mp(tmp_sorted_errors(s)+256,1);
                    end
                 else
                       fake=pred_img(tmp_sorted_err_pos(s,1),tmp_sorted_err_pos(s,2))+mp(tmp_sorted_errors(s)+256,1); 
                 end
                
                if fake<=tmplayer-1%%��ֹ���ͬʱ����λͼLB
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
LB=[LB 0 tmp_LB];%%������ѹ��
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
%�������Ϣ
%���룺
%   P����ֵ��
%   Z�����
%   length��ʵ��ʹ�ò�ֵ��
%�����
%   si����������Ϣ��������ֵ��������ֵ��-��㣬ʵ��ʹ�ò�ֵ��

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
%����ÿһ��ƽ�ƺ�Ƕ��ĵ��ӹ����ֵӳ���
%���룺
%   P����ֵ��
%   Z�����
%�����
%   mp����ֵӳ�������ֵi���Ƿ�ֵ�㣬��mp��i+256��1�� = mp��i+256��2��Ϊƽ�ƺ��ֵ,
%       ����mp��i+256��1��ΪǶ��0��Ĳ�ֵ��mp��i+256��2��ΪǶ��1��Ĳ�ֵ

mp = zeros(511,2);
P = P+256;
Z = Z+256;
for i=1:numel(P)
    if P(i) < Z(i)
        if P(i)+1 <= Z(i)-1
            mp(P(i)+1:Z(i)-1,:) = mp(P(i)+1:Z(i)-1,:) + 1;                 %ƽ��
        end
        mp(P(i),2) = mp(P(i),2) + 1;                                       %Ƕ��
    else
        if P(i)-1 >= Z(i)+1
            mp(Z(i)+1:P(i)-1,:) = mp(Z(i)+1:P(i)-1,:) - 1;                 %ƽ��
        end
        mp(P(i),2) = mp(P(i),2) - 1;                                       %Ƕ��
    end
end
for i = -255:255
    mp(i+256,:) = mp(i+256,:) + i;                                         %Ԥ���������޸������Ԥ�����
end


%%
function [pred_errors_cnn, pred_img_cnn, err_pos_cnn] = predict_one_set(img, pred_img,set,side)
%���Ԥ��ͼ���Ԥ���ֵ
%���룺
%   img��ԭͼ��
%   set��X/O����0����X��
%   side����ѡȡ�ķ�Χ��3*3��5*5��7*7��
%�����
%   pred_errors��Ԥ���ֵ����
%   pred_img��Ԥ��ͼ��
%   img_pos��Ԥ���ֵ��Ӧ��ͼ��λ��

[r,c] = size(img);
img = double(img);
pre_img = double(pred_img);
pred_img_cnn = img;
pred_img_cnn(side+1:r-side,side+1:c-side) = pre_img(side+1:r-side,side+1:c-side);
% ��ԭͼ������һ�����ϵ�����ֵ��Ԥ��ļ�����ϳ�һ�����ϵ�Ԥ��ͼ
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
%ѡ����������Ȩ�أ�����Ԥ���ֵ
%���룺
%   img��Ԥ��ֵͼ��
%   pred_errors��Ԥ���ֵ����
%   err_pos��Ԥ���ֵ��ͼ����λ��
%   set��X/O����־��0����X��
%   method,���õ�����������0�����������1����Sachnev��������2����Hwang������
%�����
%   sorted_errors��������Ԥ���ֵ����
%   sorted_err_pos��������ֵ��ͼ����λ��
%   sort_side_info�����������У�������Ż�����
%   record_info����¼���ڹ۲����Ϣ
%   cmt,��������ֵ���õ�ʱ��
%   sortsurS,������Ԥ������������������������Ӧ������ֵ
%   corr��δ�����Ԥ�����ľ���ֵ�������ֵ�����ϵ��
%   featurenumber��������

%��ʼ��
[r,c] = size(img);
img = double(img);
sorted_errors = pred_errors;
sorted_err_pos = err_pos;


if method==0
        [surS,w,cmt,sur,ms,originw]=old_feature(img,pred_errors,r,c,set,side,dm);
        %%tmp_wȨ��ֵ
        %%tmp_surδ�����Ԥ�����
        %%tmp_originwδ��������Ȩ��ֵ
        %%msֵ����������
        %%tmp_cmt���õ�ʱ��
        %%surS���еĵ�����ֵ
        sortsurS=zeros(ms+2,numel(sur));
         [~,pos] = sort(sur);
             for j = 1:numel(sur)
                 sorted_errors(1,j) = pred_errors(1,pos(j));
                 sorted_err_pos(j,1) = err_pos(pos(j),1);
                 sorted_err_pos(j,2) = err_pos(pos(j),2);
                 sortsurS(1:ms,j)=surS(1:ms,pos(j));   
                 sortsurS(ms+1,j)=sur(1,pos(j));%%�����Ķ�����ֵ
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
    error('methodѡ�����');
end

if ismember(method,[1 2 3])
%����
[~,pos] = sort(sur);
sortsurS=zeros(ms+2,numel(sur));

% ���ݸ��Ӷȴ�Сѡ��ͬԤ����������Ԥ����0��ʾ��cnnԤ������1��ʾ
[pred_errors, pred_img, err_pos, prediction_map] = prediction_select(img,pred_img_cnn,side,set,method);
record_info.prediction_map = prediction_map;

    for i = 1:numel(sur)
      sorted_errors(1,i) = pred_errors(1,pos(i));
      sorted_err_pos(i,1) = err_pos(pos(i),1);
      sorted_err_pos(i,2) = err_pos(pos(i),2);
      sortsurS(1:ms,i)=surS(1:ms,pos(i));   
      sortsurS(ms+1,i)=sur(1,pos(i));%%�����Ķ�����ֵ
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
%%����
%%img����ͼ��
%%pred_errorsԤ�����
%%r,cͼ��ĳ��Ϳ�
%%set��0��ʾx���ϣ�1��ʾo����
%%method 1��ʾSachnev3*3��ֵ����ֵ�ķ���
%%method2 ��ʾHwang3*3����
%%method3��ʾLi������
%%side��ʾѡ�������Ŀ�ȣ�ָͼ�����side���ȵ����ز�ִ��Ƕ��
t1=clock();
sur=zeros(1,numel(pred_errors));
sur1=sur;
k=1;
for row = side+1:r-side
    for col = side+1:c-side
        if mod(row+col,2) == set
            if method==1
            surround = [img(row,col-1),img(row-1,col),img(row,col+1),img(row+1,col)];%3*3��������
            dist_near = abs(diff([surround,img(row,col-1)]));%�������صĲ�ֵ
            sur1(1,k) = var(dist_near,1);
            end
            if method==2
            surround = [img(row,col-1),img(row-1,col),img(row,col+1),img(row+1,col)];%3*3��������
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
% for row = side+1:r-side                           %%ʹ���ɵ�����ֵ�����ظ�
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
%����ֵ����
sur = zeros(1,numel(pred_errors));
%12������ֵ����
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

%�ԽǵĶ��ײ�ֵ����
filterBlk = [1,0,1;0,-4,0;1,0,1];
%���ڼ����ݶȵ����ص�λ������м����ص�ƫ��
gradOffset = [-1,0,1,0;-2,1,0,1;-1,2,1,2;0,1,2,1];
gradOffset2 = [-1,-2,1,-2;-2,-1,0,-1;-1,0,1,0;0,-1,2,-1];
%�ݶȲ�ֵ
gradDiff = zeros(6,2);
gradDiff2 = zeros(6,2);
k = 1;
tt1 = clock();
for row = side+1:r-side
    for col = side+1:c-side
        if mod(row+col,2) == set
            %3*3��������
            surround = [img(row,col-1),img(row-1,col),img(row,col+1),img(row+1,col)];
            %�������صĲ�ֵ
            dist_near = abs(diff([surround,img(row,col-1)]));
                %��������������ص��ݶ�
                grad = [img(row+gradOffset(1,1),col+gradOffset(1,2))-surround(1),img(row+gradOffset(1,3),col+gradOffset(1,4))-surround(1);
                    img(row+gradOffset(2,1),col+gradOffset(2,2))-surround(2),img(row+gradOffset(2,3),col+gradOffset(2,4))-surround(2);
                    img(row+gradOffset(3,1),col+gradOffset(3,2))-surround(3),img(row+gradOffset(3,3),col+gradOffset(3,4))-surround(3);
                    img(row+gradOffset(4,1),col+gradOffset(4,2))-surround(4),img(row+gradOffset(4,3),col+gradOffset(4,4))-surround(4)];
                grad2 = [img(row+gradOffset2(1,1),col+gradOffset2(1,2))-surround(1),img(row+gradOffset2(1,3),col+gradOffset2(1,4))-surround(1);
                    img(row+gradOffset2(2,1),col+gradOffset2(2,2))-surround(2),img(row+gradOffset2(2,3),col+gradOffset2(2,4))-surround(2);
                    img(row+gradOffset2(3,1),col+gradOffset2(3,2))-surround(3),img(row+gradOffset2(3,3),col+gradOffset2(3,4))-surround(3);
                    img(row+gradOffset2(4,1),col+gradOffset2(4,2))-surround(4),img(row+gradOffset2(4,3),col+gradOffset2(4,4))-surround(4)];
                %��������������ص��ݶȲ�
                gk = 1;
                for p = 1:3
                    for q = p+1:4
                        gradDiff(gk,:) = grad(p,:) - grad(q,:);
                        gradDiff2(gk,:) = grad2(p,:) - grad2(q,:);
                        gk = gk+1;
                    end
                end
            

            %Hwang��Sachnev���õ�����ֵ
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

%����������ֵ���ŵ�[0,1000]�ķ�Χ��
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
    absErr = abs(pred_errors);%%����Ȩ��
    errNum = numel(absErr);
    lb=ones(1,10).*(-255);%%Ȩ��������
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
w = int32(w/max(abs(w))*255);%��ΧΪ-255~255
% if di
% sur=  double(originw)*surS;
% intergersur=double(w)*surS;
% else
sur = double(w)*surS;
% intergersur=double(w)*surS;
% end


k=1;
for row = side+1:r-side                           %%ʹ���ɵ�����ֵ�����ظ�
    for col = side+1:c-side
        if mod(row+col,2)==set
    sur(k)=sur(k)+(row+(col-1)*r)/10^(15);
    k=k+1;
        end
    end
end
tt2 = clock();
cmt = etime(tt2,tt1);
ms=10;%%our������ʹ�õ�����
%%
function [res] = calCorr(w,absErr,surS,k)
%����Ȩ��Ϊw������ֵΪsurSʱ��Ԥ���ֵ�����ϵ��
%���룺
%   w��Ȩ��
%   err����ֵ����
%   surS��ȫ������ֵ
%   k����ֵ����
%�����
%   res�����ϵ��
fea = double(w)*surS;
meanF = mean(fea);
res  = 0-(absErr*fea'/k - mean(absErr)*meanF)/(sqrt(mean((fea).^2) - (meanF)^2)*sqrt(mean((absErr).^2)-(mean(absErr)).^2));
%%
function [BS]  = toBitSeq(x,l)%�������x��w(i)ת��Ϊһ��lλ(��8λ)��������
%��x�õ�x����ֵ�Ķ���������
%����:
%   x������
%   l�����г���
%���:
%   BS������������
BS = zeros(1,l);
x = uint32(abs(x));
for i=1:l
    BS(1,i) = bitget(x,l-i+1);%��xת��Ϊ����������ȡ��λ��ʼ��iλ��
end
%%
function [BS] =  getWSeq(w)
%��Ȩ�ر���Ϊ01����
%���룺
%   w��Ȩ��
%   BS�����������
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




