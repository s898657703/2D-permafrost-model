function [TPP,TSS]=FMULTI(NUM_node,NNN,NN,TSS,TPP,TTT,TNN)

for m = 1:NUM_node 
    for n = 1:NUM_node   
        if (m>=n)
            II7 = m;
            II8 = n;

            
            if ((II7-II8)<=(NNN(II7+1)-NNN(II7)-1))
                TSS(II7) = TNN(NN(II7)-II7+II8);
            else
                TSS(II7) = 0;
            end
        else
            II7 = n;
            II8 = m;

            if ((II7-II8)<=(NNN(II7+1)-NNN(II7)-1))
                TSS(II7) = TNN(NN(II7)-II7+II8);
            else
                TSS(II7) = 0;
            end
        end
        if (TSS(II7)~=0)
            TPP(m) = TPP(m)+TSS(II7)*TTT(n);

        end
    end
end
end
%%