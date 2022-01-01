function [TTT]=FGAUSSIAN(NUM_node,NNN,NN,TSS,TPP,TKK)
%MULTI针对系数矩阵N,判断矩阵系数是否为0，是否要执行消去运算140页
       %GAUSSIAN  METHOD高斯迭代法124页
        for m = 1:NUM_node%m20--II
            TPP(m) = TPP(m)/TKK(NN(m));
            if (m~=NUM_node)
                for n = m+1:NUM_node%m21--JJ
                    II7 = n;
                    II8 = m;
                    
                    if ((II7-II8)<=(NNN(II7+1)-NNN(II7)-1))
                        TSS(II7) = TKK(NN(II7)-II7+II8);
                    else
                        TSS(II7) = 0;
                    end
                    if (TSS(n)~=0)
                        IAO = NN(n)-n+m;
                        TSS(n) = TKK(IAO)/TKK(NN(m));
                        for m22 = IAO+1:NN(n)%m22--IO
                            TKK(m22) = TKK(m22)-TKK(IAO)*TSS(m+m22-IAO);
                        end
                        TPP(n) = TPP(n)-TKK(IAO)*TPP(m);
                        TKK(IAO) = TSS(n);
                    end
                end
            end
            TKK(NN(m)) = 1;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for m = NUM_node-1:-1:1%m23--II
            for n = m+1:NUM_node%m24--JJ
                II7=n;
                II8=m;
                
                if ((II7-II8)<=(NNN(II7+1)-NNN(II7)-1))
                    TSS(II7) = TKK(NN(II7)-II7+II8);
                else
                    TSS(II7) = 0;
                end
                TPP(m) = TPP(m)-TSS(n)*TPP(n);
            end
        end
        
         TTT=TPP;
end
%%