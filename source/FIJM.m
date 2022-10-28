function [IJM]=FIJM(NY,NX)

for i = 1:NY-1
    for j = 1:2:2*(NX-1)-1
        if (j<2*(NX-1)-1)
            IJM((i-1)*(NX-1)*2+j,1) = NX+2+NX*(i-1)+(j-1)/2;
            IJM((i-1)*(NX-1)*2+j,2) = NX+1+NX*(i-1)+(j-1)/2;
            IJM((i-1)*(NX-1)*2+j,3) = 1+NX*(i-1)+(j-1)/2;
            IJM((i-1)*(NX-1)*2+j+1,1) = NX+2+NX*(i-1)+(j-1)/2;
            IJM((i-1)*(NX-1)*2+j+1,2) = 1+NX*(i-1)+(j-1)/2;
            IJM((i-1)*(NX-1)*2+j+1,3) = 2+NX*(i-1)+(j-1)/2;
        else
            IJM((i-1)*(NX-1)*2+j,1) = NX*2-1+NX*(i-1);
            IJM((i-1)*(NX-1)*2+j,2) = NX-1+NX*(i-1);
            IJM((i-1)*(NX-1)*2+j,3) = NX+NX*(i-1);
            IJM((i-1)*(NX-1)*2+j+1,1) = NX*2-1+NX*(i-1);
            IJM((i-1)*(NX-1)*2+j+1,2) = NX+NX*(i-1);
            IJM((i-1)*(NX-1)*2+j+1,3) = 2*NX+NX*(i-1);
        end
    end
end

end

