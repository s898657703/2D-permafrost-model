function [NUM_node,NUM_tri,NUM_upnode,NUM_bontri,NODE]=FCALCULATOR(NX,NY)
%%

NUM_node = NY*NX;% �ڵ��ܸ���
NUM_tri = (NY-1)*(NX-1)*2;%�����ε�Ԫ����
NUM_upnode = NX; %�ϱ߽����
NUM_bontri = (NX-1)+(NY-1)*2;%�߽������θ���
%%  Ϊÿ���ڵ����һ��Ψһ�ڵ��
for i = 1:NX
    for ii = 1:NY
        NODE(ii,i) = (NY-ii)*NX+i;%�ڵ��
    end
end
end