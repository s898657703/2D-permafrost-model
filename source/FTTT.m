function [TTT]=FTTT(T_T,NX,NY)
%%      ���ó�ʼ�¶�ֵ�����¶ȳ�

%���������¶ȳ��������
for i =  1:NY
    for j = 1:NX
        TTT(NX*(NY-i)+j,1) = T_T(i,j);
    end
end
end