function [NUM_node,NUM_tri,NUM_upnode,NUM_bontri,NODE]=FCALCULATOR(NX,NY)
%%

NUM_node = NY*NX;
NUM_tri = (NY-1)*(NX-1)*2;
NUM_upnode = NX;
NUM_bontri = (NX-1)+(NY-1)*2;
for i = 1:NX
    for ii = 1:NY
        NODE(ii,i) = (NY-ii)*NX+i;
    end
end
end