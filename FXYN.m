function [XYN]=FXYN(NY,NX,NUM_node,ipx,ipy)
% 为每个三角形元素分配一个唯一的元素编号。为每个元素指定三个顶点的节点号。
% IJM存储的是每个三角形元素的三个节点号，编号从左下角开始;且都是jm边位于边界

XYN = zeros(NUM_node,2);%对应 XYN(INT0,2);XYN(10000,2)
for i = 1:NY
    for j = 1:NX
        XYN(j+(i-1)*NX,2) = ipx(j);
        XYN(j+(i-1)*NX,1) = ipy(NY+1-i);
    end
end

end







