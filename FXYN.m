function [XYN]=FXYN(NY,NX,NUM_node,ipx,ipy)
% Ϊÿ��������Ԫ�ط���һ��Ψһ��Ԫ�ر�š�Ϊÿ��Ԫ��ָ����������Ľڵ�š�
% IJM�洢����ÿ��������Ԫ�ص������ڵ�ţ���Ŵ����½ǿ�ʼ;�Ҷ���jm��λ�ڱ߽�

XYN = zeros(NUM_node,2);%��Ӧ XYN(INT0,2);XYN(10000,2)
for i = 1:NY
    for j = 1:NX
        XYN(j+(i-1)*NX,2) = ipx(j);
        XYN(j+(i-1)*NX,1) = ipy(NY+1-i);
    end
end

end







