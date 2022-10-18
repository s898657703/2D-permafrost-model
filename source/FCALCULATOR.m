function [NUM_node,NUM_tri,NUM_upnode,NUM_bontri,NODE]=FCALCULATOR(NX,NY)
%%

NUM_node = NY*NX;% 节点总个数
NUM_tri = (NY-1)*(NX-1)*2;%三角形单元个数
NUM_upnode = NX; %上边界点数
NUM_bontri = (NX-1)+(NY-1)*2;%边界三角形个数
%%  为每个节点分配一个唯一节点号
for i = 1:NX
    for ii = 1:NY
        NODE(ii,i) = (NY-ii)*NX+i;%节点号
    end
end
end