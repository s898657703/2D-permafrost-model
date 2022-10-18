function [T_T]=FINI(T_1,T_2,T_3,ipx,ipy)
%%

x = [0 47.5 93.5];
T_ini1 = T_1(:,2);
y1 = T_1(:,1);

DT_y1 = 0.03;
int_ay1 = floor(y1(end));
for i = int_ay1 + 2:2:100
    y1 = [y1;i];
    T_ini1 = [T_ini1;T_ini1(end)+DT_y1*(y1(end)-y1(end-1))];
end
%%

T_ini2 = T_2(:,2);
y2 = T_2(:,1);

DT_y2 = 0.03;
int_ay2 = floor(y2(end));
for i = int_ay2 + 2:2:100
    y2 = [y2;i];
    T_ini2 = [T_ini2;T_ini2(end)+DT_y2*(y2(end)-y2(end-1))];
end
%%

T_ini3 = T_3(:,2);
y3 = T_3(:,1);

DT_y3 = 0.03;
int_ay3 = floor(y3(end));
for i = int_ay3 + 2:2:100
    y3 = [y3;i];
    T_ini3 = [T_ini3;T_ini3(end)+DT_y3*(y3(end)-y3(end-1))];
end
%%
T_init(:,1)=interp1(y1,T_ini1,ipy,'Linear');%
T_init(:,2)=interp1(y2,T_ini2,ipy,'Linear');%
T_init(:,3)=interp1(y3,T_ini3,ipy,'Linear');%

%%
T_T = interp2(x,ipy,T_init,ipx',ipy,'Linear');

end
