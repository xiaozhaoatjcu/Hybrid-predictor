function[PZ_infor]=simply_PZ_infor(Map_PZ,PZ_infor,Map_layer_infor)
%%ÿһ����ֱ��ͼ����ʹ���˶��ٶԷ�ֵ�����ֵ�����Ϣ����ѹ��
%%������ֱ��ͼ��ѡ���˴��ڵ���������ֵ�����ֵ�㣬���ٽ��п���ѹ����Ϣͳһʹ��80λ���б�ʾ
    key=sum(Map_layer_infor>2);
    a=numel(Map_layer_infor);
    if key==0
    else
        PZ_infor=[PZ_infor 0];%%��ѹ��
        for i=1:a
        PZ_layer=Map_PZ(i,end);
        PZ_layer=toBitSeq(PZ_layer,5);%%ÿ��ֱ��ͼʹ���˶��ٸ���ֵ�㣬��λ�ɴ��31λ
        PZ_infor=[PZ_infor PZ_layer];
        end
        return
    end
    key=sum(Map_layer_infor==1);
    if key==0
            %%ѹ����2+16λ
        PZ_infor=[PZ_infor 1 0];%%����ѹ����ֻ��0��2��������
        [~,pos]=find(Map_layer_infor==2);
        tmpdata=zeros(1,a);
            for i=1:numel(pos)
                tmpdata(pos(i))=1;
            end
       PZ_infor=[PZ_infor tmpdata];
    else
       %%ѹ����2+32λ
       PZ_infor=[PZ_infor 1 1];%%����ѹ������0 1 2��������
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