function [TPP,TSS]=FMULTI(NUM_node,NNN,NN,TSS,TPP,TTT,TNN)
%MULTI针对系数矩阵N,判断矩阵系数是否为0，是否要执行消去运算140页
for m = 1:NUM_node %m4--L3
    for n = 1:NUM_node   %n4--L4
        if (m>=n)
            II7 = m;
            II8 = n;
            %JUDGE判断是否执行消去计算
            
            if ((II7-II8)<=(NNN(II7+1)-NNN(II7)-1))
                TSS(II7) = TNN(NN(II7)-II7+II8);
            else
                TSS(II7) = 0;
            end
        else
            II7 = n;
            II8 = m;
            %JUDGEJUDGE判断是否执行消去计算
            if ((II7-II8)<=(NNN(II7+1)-NNN(II7)-1))
                TSS(II7) = TNN(NN(II7)-II7+II8);
            else
                TSS(II7) = 0;
            end
        end
        if (TSS(II7)~=0)
            TPP(m) = TPP(m)+TSS(II7)*TTT(n);
            %采用C-N格式计算的基本方程右侧的一个乘法矩阵计算方法
        end
    end
end
end
%%