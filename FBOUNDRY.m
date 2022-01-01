function [NODE_up,Q2]=FBOUNDRY(T_T,ipx,ipy,GST,day_run,k)
%ȷ���߽�������Ԫ�أ����ǵ����������ͣ����Ա߽�ֵ����������
%ģ����������е����ݶȸ���ʵ��ֵ������
TIME_s=3600.0 * 24;

[~,nsf] = find(ipx==47.5);%�ҵ�47.5��SFGT���λ������λ��
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


%%      �ϱ߽�Ľڵ��
for i = 1:NX%%%�����ϱ߽��
    NODE_up(i,1) = (NY-1)*NX+i;
    Tup_gradient=0.06;%���ݶ�
    %      Tup_gradient=(T_T(end,i)-T_T(end-1,i))/(ipy(end,1)-ipy(end-1,1));%���ݶ�
    heatflux_lower(1,i)=-Tup_gradient*k(end,i)*TIME_s;%������
    
end
%%
%�±߽�ı߽�������Ԫ�ر�ź�����
for i = 1:NX-2
    Q2(2*i,1) = heatflux_lower(1,i);
end

%�±߽�����Ҳ�ı߽�������Ԫ�ر��
Q2(2*(NX-1)-1,1) = heatflux_lower(1,NX-1);

%��߽�������Ԫ�ر�ź�����
for i = 1:NY-1
    Tleft_gradient=(T_T(i,1)-T_T(i,nsf))/47.5;%���ݶ�
    heatflux_left=-Tleft_gradient*k(i,1)*TIME_s;%������
    Q2(1+(i-1)*2*(NX-1)) =heatflux_left;
end

%�ұ߽�����Ҳ�ı߽�������Ԫ�ر�ź�����
for i = 1:NY-1
    Tright_gradient=(T_T(i,end)-T_T(i,nsf))/47.5;%���ݶ�
    heatflux_right=-Tright_gradient*k(i,end)*TIME_s;%������
    Q2(2*(NX-1)*i) = heatflux_right;%
end

end