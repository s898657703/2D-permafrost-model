function [NX,NY,T_T]=FINI(T_1,T_3,T_2,ipx,ipy)
%%


%%    
T_ini1 = T_1(:,2);
ay1 = T_1(:,1);
DT_y1 = (T_ini1(end,1)-T_ini1(end-1,1))/(ay1(end,1)-ay1(end-1,1)); 
int_ay1 = floor(ay1(end));
for i = int_ay1 + 2:2:100
    ay1 = [ay1;i];
    T_ini1 = [T_ini1;T_ini1(end)+DT_y1*(ay1(end)-ay1(end-1))];
end

T_ini2 = T_2(:,2);
ay2 = T_2(:,1);
  DT_y2 = (T_ini2(end,1)-T_ini2(end-1,1))/(ay2(end,1)-ay2(end-1,1));

int_ay2 = floor(ay2(end));
for i = int_ay2 + 2:2:100
    ay2 = [ay2;i];
    T_ini2 = [T_ini2;T_ini2(end)+DT_y2*(ay2(end)-ay2(end-1))];
end
%%

T_ini3 = T_3(:,2);
ay3 = T_3(:,1);
DT_y3 = (T_ini3(end,1)-T_ini3(end-1,1))/(ay3(end,1)-ay3(end-1,1)); 
int_ay3 = floor(ay3(end));
for i = int_ay3 + 2:2:100
    ay3 = [ay3;i];
    T_ini3 = [T_ini3;T_ini3(end)+DT_y3*(ay3(end)-ay3(end-1))];
end
%%
T_init(:,1)=interp1(ay1,T_ini1,ipy,'Linear');%
T_init(:,2)=interp1(ay2,T_ini2,ipy,'Linear');%
T_init(:,3)=interp1(ay3,T_ini3,ipy,'Linear');%
%%

[~,	NX] = size(ipx); 
[NY,~] = size(ipy);
[~,nrr] = find(ipx==47.5);
for i = 1:NY
    DT_x1 = (T_init(i,2)-T_init(i,1))/47.5;
    DT_x2 = (T_init(i,3)-T_init(i,2))/46;
    T_T(i,1)=T_init(i,1);
    for ii=2:nrr-1
        T_T(i,ii)=T_init(i,1)+  DT_x1 *(ipx(1,ii)-ipx(1,1));
    end
     T_T(i,nrr)=T_init(i,2);
    for ii=nrr+1:NX-1
        T_T(i,ii)=T_init(i,2)+  DT_x2 *(ipx(1,ii)-ipx(1,nrr));
    end
    T_T(i,NX)=T_init(i,3);
end

end