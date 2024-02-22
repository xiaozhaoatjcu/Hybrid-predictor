function [pred_errors, pred_img, err_pos,prediction_map] = prediction_select(img,pred_img_cnn,side,set,method)
%%根据复杂度大小设置阈值，对于不同的图像块选择不同预测器进行预测，得到的最终预测图
%%输入
%   img：原图像
%   pred_img_cnn：利用cnn预测得出的预测图
%   set：X/O集标志，0代表X集
%   method：计算复杂度不同方法
%%输出
%   pred_img：根据复杂度大小选择不同预测器进行预测得到的最终预测图
%   pred_errors：预测差值数组
%   err_pos：预测差值对应的图像位置
%   prediction_map：记录每一个图像块所使用的预测器

[r,c] = size(img);
img = double(img);
pred_img_cnn = double(pred_img_cnn);% 利用cnn预测得到的预测图
pred_img_diamond = double(img);% 利用菱形预测得到的预测图
pred_errors_diamond = zeros(1,(r-side*2)*(c-side*2)/2);% 利用菱形预测得到的预测误差
err_pos_diamond = zeros(numel(pred_errors_diamond),2);% 利用菱形预测得到的预测误差对应的图像位置

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

pred_errors_cnn = zeros(1,(r-side*2)*(c-side*2)/2);% 利用cnn预测得到的预测误差
err_pos_cnn = zeros(numel(pred_errors_cnn),2);% 利用cnn预测得到的预测误差对应的图像位置

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

% 对原图像计算复杂度
sur_temp=zeros(r,c);  % 计算得出的复杂度值放在与原图像相同的位置
% k=1;
for row = side+1:r-side
    for col = side+1:c-side
        if mod(row+col,2) == set
            if method==1
            surround = [img(row,col-1),img(row-1,col),img(row,col+1),img(row+1,col)];%3*3邻域像素
            dist_near = abs(diff([surround,img(row,col-1)]));%邻域像素的差值
            sur_temp(row,col) = var(dist_near,1);
            end
            if method==2
            surround = [img(row,col-1),img(row-1,col),img(row,col+1),img(row+1,col)];%3*3邻域像素
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

% 阈值的选取，根据复杂度不同而变化，每张图的阈值选取都不同
s_max = max(max(sur_temp));
if s_max>=700 && s_max<800
    s = round(s_max*0.157123);  % 0.151955，0.157123
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

prediction_map = zeros(1,(r-side*2)*(c-side*2)/2); % 记录每一个图像块所使用的预测器
pred_errors = zeros(1,(r-side*2)*(c-side*2)/2); % 最终预测差值序列
err_pos = zeros(numel(pred_errors_cnn),2); % 预测差值对应的图像位置
k = 1;
for row = side+1:r-side
    for col = side+1:c-side
        if mod(row+col,2) == set
            if sur_temp(row,col)>s
                pred_img(row,col) = pred_img_cnn(row,col); % cnn预测器用1表示
                prediction_map(k) = 1;
                pred_errors(k) = img(row,col) - pred_img(row,col);
                err_pos(k,1) = row;
                err_pos(k,2) = col;
                k = k + 1;
            else
                pred_img(row,col) = pred_img_diamond(row,col); % 菱形预测用0表示
                prediction_map(k) = 0;
                pred_errors(k) = img(row,col) - pred_img(row,col);
                err_pos(k,1) = row;
                err_pos(k,2) = col;
                k = k + 1;
            end
        end
    end
end

                
                
                
            








