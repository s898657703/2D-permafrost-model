function [NODE_up,Q2]=FBOUNDRY(ipx,GST_day,TIME_s,wide_left,wide_right,NX,NY,heat_flux)
%Determine the number and type of boundary triangle elements and set corresponding boundary conditions

n1 =  find(ipx>=wide_left,1,'first');
n2 =  find(ipx>=wide_right,1,'first');

for ii=1:n1
    GST1(1,ii) =  GST_day(1,1);
end
for ii=n1+1:n2
    GST1(1,ii) =  GST_day(1,2);
end

for ii=n2+1:NX
    GST1(1,ii) =  GST_day(1,3);
end



% upper boundary
NODE_up(:,2)=GST1';
for i = 1:NX
    NODE_up(i,1) = (NY-1)*NX+i;
end
% lower boundary
for i = 1:NX-2
    Q2(2*i,1) = -heat_flux*TIME_s;
end
Q2(2*(NX-1)-1,1) = -heat_flux*TIME_s;
% left boundary
for i = 1:NY-1

    Q2(1+(i-1)*2*(NX-1)) =0;
end
% right boundary
for i = 1:NY-1
    Q2(2*(NX-1)*i) = 0;
end
end