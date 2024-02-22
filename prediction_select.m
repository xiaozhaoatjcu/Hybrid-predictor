function [pred_errors, pred_img, err_pos,prediction_map] = prediction_select(img,pred_img_cnn,side,set,method)
%%���ݸ��Ӷȴ�С������ֵ�����ڲ�ͬ��ͼ���ѡ��ͬԤ��������Ԥ�⣬�õ�������Ԥ��ͼ
%%����
%   img��ԭͼ��
%   pred_img_cnn������cnnԤ��ó���Ԥ��ͼ
%   set��X/O����־��0����X��
%   method�����㸴�ӶȲ�ͬ����
%%���
%   pred_img�����ݸ��Ӷȴ�Сѡ��ͬԤ��������Ԥ��õ�������Ԥ��ͼ
%   pred_errors��Ԥ���ֵ����
%   err_pos��Ԥ���ֵ��Ӧ��ͼ��λ��
%   prediction_map����¼ÿһ��ͼ�����ʹ�õ�Ԥ����

[r,c] = size(img);
img = double(img);
pred_img_cnn = double(pred_img_cnn);% ����cnnԤ��õ���Ԥ��ͼ
pred_img_diamond = double(img);% ��������Ԥ��õ���Ԥ��ͼ
pred_errors_diamond = zeros(1,(r-side*2)*(c-side*2)/2);% ��������Ԥ��õ���Ԥ�����
err_pos_diamond = zeros(numel(pred_errors_diamond),2);% ��������Ԥ��õ���Ԥ������Ӧ��ͼ��λ��

k = 1;
for row = side+1:r-side
    for col = side+1:c-side
        if mod(row + col, 2) == set
            pred_img_diamond(row,col) = round(mean([img(row-1,col),img(row,col+1),img(row+1,col),img(row,col-1)]));
            pred_errors_diamond(k) = img(row,col) - pred_img_diamond(row,col);
            err_pos_diamond(k,1) = row;
            err_pos_diamond(k,2) = col;
            k = k + 1;
        end
    end
end

pred_errors_cnn = zeros(1,(r-side*2)*(c-side*2)/2);% ����cnnԤ��õ���Ԥ�����
err_pos_cnn = zeros(numel(pred_errors_cnn),2);% ����cnnԤ��õ���Ԥ������Ӧ��ͼ��λ��

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

% ��ԭͼ����㸴�Ӷ�
sur_temp=zeros(r,c);  % ����ó��ĸ��Ӷ�ֵ������ԭͼ����ͬ��λ��
% k=1;
for row = side+1:r-side
    for col = side+1:c-side
        if mod(row+col,2) == set
            if method==1
            surround = [img(row,col-1),img(row-1,col),img(row,col+1),img(row+1,col)];%3*3��������
            dist_near = abs(diff([surround,img(row,col-1)]));%�������صĲ�ֵ
            sur_temp(row,col) = var(dist_near,1);
            end
            if method==2
            surround = [img(row,col-1),img(row-1,col),img(row,col+1),img(row+1,col)];%3*3��������
            sur_temp(row,col)= var(surround,1);
            end
            if method==3
                surround3=[abs(img(row,col-1)-img(row+1,col-1)),abs(img(row+1,col-1)-img(row+2,col-1)),abs(img(row+1,col)-img(row+2,col))...
             ,abs(img(row,col+1)-img(row+1,col+1)),abs(img(row+1,col+1)-img(row+2,col+1)),abs(img(row-1,col+2)-img(row,col+2))...
             ,abs(img(row,col+2)-img(row+1,col+2)),abs(img(row+1,col+2)-img(row+2,col+2)),abs(img(row,col+1)-img(row,col+2))...
             ,abs(img(row+1,col-1)-img(row+1,col)),abs(img(row+1,col)-img(row+1,col+1)),abs(img(row+1,col+1)-img(row+1,col+2))...
             ,abs(img(row+2,col-1)-img(row+2,col)),abs(img(row+2,col)-img(row+2,col+1)),abs(img(row+2,col+1)-img(row+2,col+2))];
            sur_temp(row,col)=sum(surround3);
            end
%             k=k+1;
        end
    end
end

pred_img = double(img);

% ��ֵ��ѡȡ�����ݸ��ӶȲ�ͬ���仯��ÿ��ͼ����ֵѡȡ����ͬ
s_max = max(max(sur_temp));
if s_max>=700 && s_max<800
    s = round(s_max*0.157123);  % 0.151955��0.157123
elseif s_max>=800 && s_max<900
    s = round(s_max*0.06536);
elseif s_max>=900 && s_max<950
    s = round(s_max*0.105);
elseif s_max>=950 && s_max<1000
    s = round(s_max*0.053625);
elseif s_max>=1000
    s = round(s_max*0.1567);
else 
    s = s_max;
end   

prediction_map = zeros(1,(r-side*2)*(c-side*2)/2); % ��¼ÿһ��ͼ�����ʹ�õ�Ԥ����
pred_errors = zeros(1,(r-side*2)*(c-side*2)/2); % ����Ԥ���ֵ����
err_pos = zeros(numel(pred_errors_cnn),2); % Ԥ���ֵ��Ӧ��ͼ��λ��
k = 1;
for row = side+1:r-side
    for col = side+1:c-side
        if mod(row+col,2) == set
            if sur_temp(row,col)>s
                pred_img(row,col) = pred_img_cnn(row,col); % cnnԤ������1��ʾ
                prediction_map(k) = 1;
                pred_errors(k) = img(row,col) - pred_img(row,col);
                err_pos(k,1) = row;
                err_pos(k,2) = col;
                k = k + 1;
            else
                pred_img(row,col) = pred_img_diamond(row,col); % ����Ԥ����0��ʾ
                prediction_map(k) = 0;
                pred_errors(k) = img(row,col) - pred_img(row,col);
                err_pos(k,1) = row;
                err_pos(k,2) = col;
                k = k + 1;
            end
        end
    end
end

                
                
                
            








