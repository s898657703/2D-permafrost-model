function [TKK,TNN,TPP]=FCRANK(NUM_node,TPP,NUM_max,TKK,TNN)
                 
        for i = 1:NUM_node
            TPP(i) = 2*TPP(i);
        end
        TS = 1;
        TEE = zeros(NUM_max,1);
        for i = 1:NUM_max  
            TEE(i) = TKK(i);
            TKK(i) = TKK(i)+TNN(i)*2/TS;
            TNN(i) = TNN(i)*2/TS-TEE(i);
        end

end