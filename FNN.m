function [NN,NUM_max,NNN]=FNN(IJM,NUM_node,NUM_tri)
%地址子程序137页；确定矢量存储器中系数矩阵的每个主要对角元素的地址和最大存储数量
NN = zeros(NUM_node,1);
NNN = zeros(NUM_node+1,1);
for i = 1:NUM_tri%!求K系数矩阵中非0元素宽度
    IJME = min(IJM(i,:));
    if ( (NN(IJM(i,1),1))<=(IJM(i,1)-IJME))
        NN(IJM(i,1),1)=IJM(i,1)-IJME+1;
    end
    if ((NN(IJM(i,2),1))<=(IJM(i,2)-IJME))
        NN(IJM(i,2),1)=IJM(i,2)-IJME+1;
    end
    if ((NN(IJM(i,3),1))<=(IJM(i,3)-IJME))
        NN(IJM(i,3),1)=IJM(i,3)-IJME+1;
    end
end
%NN存储的是非0元素的宽度，经过下式计算后NN存储的是每个节点系数矩阵的地址
for j = 2:1:NUM_node
    NN(j,1) = NN(j,1) + NN(j-1,1);
end
NUM_max = NN(NUM_node);%最大存储数

for i = 2:NUM_node+1
NNN(i)=NN(i-1);
end







