function [TKK,TNN,TPP]=FCRANK(NUM_node,TPP,NUM_max,TKK,TNN)
                  % %CRANK-NICOLSON
        for i = 1:NUM_node%节点数
            TPP(i) = 2*TPP(i);
        end
        TS = 1;%时间步长
        TEE = zeros(NUM_max,1);
        for i = 1:NUM_max  %根据C-N方式生成新的计算矩阵
            TEE(i) = TKK(i);
            TKK(i) = TKK(i)+TNN(i)*2/TS;
            TNN(i) = TNN(i)*2/TS-TEE(i);
        end

end