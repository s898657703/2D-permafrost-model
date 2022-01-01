function [ipx,ipy]=Fipxy

ipx0=(0:0.2:1);
ipx1=(1.5:0.5:2);
ipx2=(4:2:46);
ipx3=(46.6:0.2:47.4);
ipx4=47.5;
ipx5=(47.6:0.2:48);
ipx6=(48.5:0.5:49);
ipx7=(50:2:92);
ipx8=(92.6:0.2:93.4);
ipx9=93.5;
ipx=[ipx0 ipx1 ipx2 ipx3 ipx4 ipx5 ipx6 ipx7 ipx8 ipx9];
%%
ipy1=(0.1:0.1:15)';
ipy2=(15.5:0.5:22)';
ipy3=(24:2:30)';
ipy4=(35:5:100)';
ipy=[ipy1;ipy2;ipy3;ipy4];
end