function[PZ_infor]=simply_PZ_infor(Map_PZ,PZ_infor,Map_layer_infor)
%%每一个子直方图各自使用了多少对峰值点和零值点的信息进行压缩
%%对于子直方图中选用了大于等于三个峰值点和零值点，不再进行考虑压缩信息统一使用80位进行表示
    key=sum(Map_layer_infor>2);
    a=numel(Map_layer_infor);
    if key==0
    else
        PZ_infor=[PZ_infor 0];%%不压缩
        for i=1:a
        PZ_layer=Map_PZ(i,end);
        PZ_layer=toBitSeq(PZ_layer,5);%%每个直方图使用了多少个峰值点，五位可存放31位
        PZ_infor=[PZ_infor PZ_layer];
        end
        return
    end
    key=sum(Map_layer_infor==1);
    if key==0
            %%压缩成2+16位
        PZ_infor=[PZ_infor 1 0];%%进行压缩且只有0和2两种数字
        [~,pos]=find(Map_layer_infor==2);
        tmpdata=zeros(1,a);
            for i=1:numel(pos)
                tmpdata(pos(i))=1;
            end
       PZ_infor=[PZ_infor tmpdata];
    else
       %%压缩成2+32位
       PZ_infor=[PZ_infor 1 1];%%进行压缩且有0 1 2三种数字
       [~,pos]=find(Map_layer_infor==2);
        tmpdata=zeros(1,a);
            for i=1:numel(pos)
                tmpdata(pos(i))=1;
            end
       PZ_infor=[PZ_infor tmpdata];
       [~,pos]=find(Map_layer_infor==1);
        tmpdata=zeros(1,a);
            for i=1:numel(pos)
                tmpdata(pos(i))=1;
            end
       PZ_infor=[PZ_infor tmpdata];
    end
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