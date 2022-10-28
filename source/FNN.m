function [NN,NUM_max,NNN]=FNN(IJM,NUM_node,NUM_tri)

NN = zeros(NUM_node,1);
NNN = zeros(NUM_node+1,1);
for i = 1:NUM_tri
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

for j = 2:1:NUM_node
    NN(j,1) = NN(j,1) + NN(j-1,1);
end
NUM_max = NN(NUM_node);

for i = 2:NUM_node+1
NNN(i)=NN(i-1);
end







