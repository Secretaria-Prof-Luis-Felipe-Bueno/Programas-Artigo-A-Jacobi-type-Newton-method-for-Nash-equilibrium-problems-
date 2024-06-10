rtx = linspace(-1.9,1.9,10);
rty = linspace(-1.9,1.9,10);

restrings = cell(10,10);
valfs1 = zeros(10,10);
valfs2 = zeros(10,10);
valfsjacobi1 = zeros(10,10); restt=zeros(10,10);
valfsjacobi2 = zeros(10,10); restiter=zeros(10,10);
meuserros = zeros(10,10);
meuserrosjacobi = zeros(10,10);
contmenor=0; contmenorjacobi=0;
%errosyuan= zeros(10,10);
for outi=1:10
    for outj=1:10
        %x0 = [tx(i);ty(j);tx(11-i);ty(11-j)];
        x0 = [-2+4*rand(1,1);-2+4*rand(1,1);-2+4*rand(1,1);-2+4*rand(1,1)];
        Geral
        %restringsjacobi{outi,outj}=strcat(num2str(xjacobi(1)),'; ',num2str(xjacobi(2)),'; ',num2str(xjacobi(3)),'; ',num2str(xjacobi(4)),'; ');
        restringsnewton{outi,outj}=strcat(num2str(x0(1)),'; ',num2str(x0(2)),'; ',num2str(x0(3)),'; ',num2str(x0(4)),'; ');
        restiter(outi,outj)=k;
        restt(outi,outj)=maint;
        %restrings{outi,outj}=strcat(num2str(x0(1)),'; ',num2str(x0(2)));
        %valfs1(outi,outj) = min(eigs(Feval(x0,A1,b1,c1,2,'aplbuenovdd',1,r1,ro1)));
        %valfs2(outi,outj) = min(eigs(Feval(x0,A2,b2,c2,2,'aplbuenovdd',2,r2,ro2)));
        %valfsjacobi1(outi,outj) = min(eigs(Feval(xjacobi,A1,b1,c1,2,'aplbuenovdd',1,r1,ro1)));
        %valfsjacobi2(outi,outj) = min(eigs(Feval(xjacobi,A2,b2,c2,2,'aplbuenovdd',2,r2,ro2)));
        %if (valfs1(outi,outj)<1e-6)||(valfs2(outi,outj)<1e-6)
        %   contmenor =contmenor+1;
        %end
        %if (valfsjacobi1(outi,outj)<1e-6)||(valfsjacobi2(outi,outj)<1e-6)
        %   contmenorjacobi =contmenorjacobi+1;
        %end
        %meuserros(outi,outj) = erro;
        %meuserrosjacobi(outi,outj) = ejacobi;
        %xyuan0 = x0+1e-4*rand(4,1);
        %[ retx,retk,t,erroyuan ] = Yuan_Varios_Fun(xyuan0,n,m,A1,b1,c1,A2,b2,c2,mystr,r1,r2,ro1,ro2);
        %restrings{i,j}=strcat(num2str(retx(1)),'; ',num2str(retx(2)),'; ',num2str(retx(3)),'; ',num2str(retx(4)),'; ');
        %errosyuan(i,j) = erroyuan;
        clear x0;
    end
end

        