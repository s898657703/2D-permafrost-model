function [NODE_up,Q2]=FBOUNDRY(T_T,ipx,ipy,GST,day_run,k)
%确定边界三角形元素，他们的数量和类型，并对边界值设置热流量
%模拟参数过程中地温梯度根据实测值计算获得
TIME_s=3600.0 * 24;

[~,nsf] = find(ipx==47.5);%找到47.5即SFGT钻孔位置所在位置
[~,	NX] = size(ipx);
[NY,~] = size(ipy);
[~,nrr] = find(ipx==47.5);
GST(:,2)=GST(:,2)+0.5;
DT_1 = (GST(day_run,2)-GST(day_run,1))/47.5;
DT_2 = (GST(day_run,3)-GST(day_run,2))/46;

NODE_up(1,2)=GST(day_run,1);

for ii=2:nrr-1
    
    NODE_up(ii,2)=GST(day_run,1)+  DT_1 *(ipx(1,ii)-ipx(1,1));
    
end
NODE_up(nrr,2)=GST(day_run,2);

for ii=nrr+1:NX-1
    NODE_up(ii,2)=GST(day_run,2)+  DT_2 *(ipx(1,ii)-ipx(1,nrr));
    
end
NODE_up(NX,2)=GST(day_run,3);


%%      上边界的节点号
for i = 1:NX%%%定义上边界点
    NODE_up(i,1) = (NY-1)*NX+i;
    Tup_gradient=0.06;%热梯度
    %      Tup_gradient=(T_T(end,i)-T_T(end-1,i))/(ipy(end,1)-ipy(end-1,1));%热梯度
    heatflux_lower(1,i)=-Tup_gradient*k(end,i)*TIME_s;%热流量
    
end
%%
%下边界的边界三角形元素编号和数量
for i = 1:NX-2
    Q2(2*i,1) = heatflux_lower(1,i);
end

%下边界的最右侧的边界三角形元素编号
Q2(2*(NX-1)-1,1) = heatflux_lower(1,NX-1);

%左边界三角形元素编号和数量
for i = 1:NY-1
    Tleft_gradient=(T_T(i,1)-T_T(i,nsf))/47.5;%热梯度
    heatflux_left=-Tleft_gradient*k(i,1)*TIME_s;%热流量
    Q2(1+(i-1)*2*(NX-1)) =heatflux_left;
end

%右边界的最右侧的边界三角形元素编号和数量
for i = 1:NY-1
    Tright_gradient=(T_T(i,end)-T_T(i,nsf))/47.5;%热梯度
    heatflux_right=-Tright_gradient*k(i,end)*TIME_s;%热流量
    Q2(2*(NX-1)*i) = heatflux_right;%
end

end