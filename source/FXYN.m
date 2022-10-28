function [XYN]=FXYN(NY,NX,NUM_node,ipx,ipy)

XYN = zeros(NUM_node,2);
for i = 1:NY
    for j = 1:NX
        XYN(j+(i-1)*NX,2) = ipx(j);
        XYN(j+(i-1)*NX,1) = ipy(NY+1-i);
    end
end

end







