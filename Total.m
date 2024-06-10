%ordem:
% testepaper
% problemayuan
% cubico
% quartico
% quadratica_yuan
% quadtestes1
clear

id =  'MATLAB:eigs:SigmaChangedToSA';
warning('off',id)



%%%%%%%quartico %%%%%%%


alfa = 1/2; 
tol = 1e-4; kmax=1000;
minvaltol = 1;
maxvaltol=1e10;
%%30 problemas de n=m=10
n=1; m=1;
mystr = 'quartico';
x0=ones(n+m,1);
kvec=zeros(30,1); %guardar iter
tvec=zeros(30,1); %guardar t
evec=zeros(30,1); %guardar erro
xvec=zeros(30,2); %guardar ponto (soh qd n=m=1)

kvecjacobi=zeros(30,1); kvecintern=kvec; kvecdois=kvec; kvecyuan = kvec;
tvecjacobi=zeros(30,1); tvecintern=tvec; tvecdois=tvec; tvecyuan = tvec;
evecjacobi=zeros(30,1); evecintern=evec; evecdois=evec; evecyuan = evec;
xvecjacobi=zeros(30,2); xvecintern=xvec; xvecdois=xvec; xvecyuan = xvec;

A1 = [];
A2 = [];
                
for i=1:30
    nustr = num2str(i);   
    b1 = (-1.5+i/10)*ones(n,1); b2=b1;
    c1 = (1.5-i/10)*ones(m,1); c2=c1;
    flaghess = 'ifdois'; multhess = 1;
    Geral
    kvec(i) = k;
    tvec(i) = maint;
    evec(i) = erro;
    xvec(i,1) = x0(1);
    xvec(i,2) = x0(2);
    
    %kvecjacobi(i) = kjacobi;
    %tvecjacobi(i) = tjacobi;
    %evecjacobi(i) = ejacobi;
    %xvecjacobi(i,1) = xjacobi(1);
    %xvecjacobi(i,2) = xjacobi(2);
    
    flaghess = 'jacobiintern'; x0=ones(n+m,1);
    Geral
    kvecintern(i) = k;
    tvecintern(i) = maint;
    evecintern(i) = erro;
    xvecintern(i,1) = x0(1);
    xvecintern(i,2) = x0(2);
    
    flaghess = 'jacobidois'; x0=ones(n+m,1);
    Geral
    kvecdois(i) = k;
    tvecdois(i) = maint;
    evecdois(i) = erro;
    xvecdois(i,1) = x0(1);
    xvecdois(i,2) = x0(2);
    
    x0=ones(n+m,1);
    [xvecret,kvecyuan(i),tvecyuan(i),evecyuan(i)] = Yuan_Varios_Fun(x0,n,m,A1,b1,c1,A2,b2,c2,mystr);
    xvecyuan(i,1)=xvecret(1);
    xvecyuan(i,2)=xvecret(2);
      
    asalvarstr = strcat('.\Salvar\resultsquartico',num2str(n),num2str(m),nustr);
    save(asalvarstr,'kvec','tvec','evec','xvec','kvecjacobi','tvecjacobi','evecjacobi','xvecjacobi','kvecintern','tvecintern','evecintern','xvecintern','kvecdois','tvecdois','evecdois','xvecdois','kvecyuan','tvecyuan','evecyuan','xvecyuan');
end

clear

erroaqui


%%%%%% quadratica_yuan %%%%%%


alfa = 1/2; 
tol = 1e-4; kmax=1000;
minvaltol = 1;
maxvaltol=1e10;
%30 problemas de n=m=10
n=1; m=1;
mystr = 'quadratica_yuan';
x0=ones(n+m,1);
kvec=zeros(30,1); %guardar iter
tvec=zeros(30,1); %guardar t
evec=zeros(30,1); %guardar erro
xvec=zeros(30,2); %guardar ponto (soh qd n=m=1)

kvecjacobi=zeros(30,1); kvecintern=kvec; kvecdois=kvec; kvecyuan = kvec;
tvecjacobi=zeros(30,1); tvecintern=tvec; tvecdois=tvec; tvecyuan = tvec;
evecjacobi=zeros(30,1); evecintern=evec; evecdois=evec; evecyuan = evec;
xvecjacobi=zeros(30,2); xvecintern=xvec; xvecdois=xvec; xvecyuan = xvec;

A1 = [];
A2 = [];
                
for i=1:30
    nustr = num2str(i);   
    b1 = (-1.5+i/10)*ones(n,1); b2=b1;
    c1 = (1.5-i/10)*ones(m,1); c2=c1;
    flaghess = 'ifdois'; multhess = 1;
    Geral
    kvec(i) = k;
    tvec(i) = maint;
    evec(i) = erro;
    xvec(i,1) = x0(1);
    xvec(i,2) = x0(2);
    
    kvecjacobi(i) = kjacobi;
    tvecjacobi(i) = tjacobi;
    evecjacobi(i) = ejacobi;
    xvecjacobi(i,1) = xjacobi(1);
    xvecjacobi(i,2) = xjacobi(2);
    
    flaghess = 'jacobiintern'; x0=ones(n+m,1);
    Geral
    kvecintern(i) = k;
    tvecintern(i) = maint;
    evecintern(i) = erro;
    xvecintern(i,1) = x0(1);
    xvecintern(i,2) = x0(2);
    
    flaghess = 'jacobidois'; x0=ones(n+m,1);
    Geral
    kvecdois(i) = k;
    tvecdois(i) = maint;
    evecdois(i) = erro;
    xvecdois(i,1) = x0(1);
    xvecdois(i,2) = x0(2);
    
    x0=ones(n+m,1);
    [xvecret,kvecyuan(i),tvecyuan(i),evecyuan(i)] = Yuan_Varios_Fun(x0,n,m,A1,b1,c1,A2,b2,c2,mystr);
    xvecyuan(i,1)=xvecret(1);
    xvecyuan(i,2)=xvecret(2);
      
    asalvarstr = strcat('.\Salvar\resultsquadratica_yuan',num2str(n),num2str(m),nustr);
    save(asalvarstr,'kvec','tvec','evec','xvec','kvecjacobi','tvecjacobi','evecjacobi','xvecjacobi','kvecintern','tvecintern','evecintern','xvecintern','kvecdois','tvecdois','evecdois','xvecdois','kvecyuan','tvecyuan','evecyuan','xvecyuan');
end


clear
erroaqui

%%%%%% quadratica yuan %%%%%%%

alfa = 1/2; 
tol = 1e-4; kmax=1000;
minvaltol = 1;
maxvaltol=1e10;
%30 problemas de n=m=10
n=10; m=10;
mystr = 'quadratica_yuan';
x0=ones(n+m,1);
kvec=zeros(30,1); %guardar iter
tvec=zeros(30,1); %guardar t
evec=zeros(30,1); %guardar erro

kvecjacobi=zeros(30,1); kvecintern=kvec; kvecdois=kvec; kvecyuan = kvec;
tvecjacobi=zeros(30,1); tvecintern=tvec; tvecdois=tvec; tvecyuan = tvec;
evecjacobi=zeros(30,1); evecintern=evec; evecdois=evec; evecyuan = evec;

A1 = [];
A2 = [];
                
for i=1:30
    nustr = num2str(i);   
    b1 = (-1.5+i/10)*ones(n,1); b2=b1;
    c1 = (1.5-i/10)*ones(m,1); c2=c1;
    flaghess = 'ifdois'; multhess = 1;
    Geral
    kvec(i) = k;
    tvec(i) = maint;
    evec(i) = erro;
    
    kvecjacobi(i) = kjacobi;
    tvecjacobi(i) = tjacobi;
    evecjacobi(i) = ejacobi;
    
    flaghess = 'jacobiintern'; x0=ones(n+m,1);
    Geral
    kvecintern(i) = k;
    tvecintern(i) = maint;
    evecintern(i) = erro;
    
    flaghess = 'jacobidois'; x0=ones(n+m,1);
    Geral
    kvecdois(i) = k;
    tvecdois(i) = maint;
    evecdois(i) = erro;
    
    x0=ones(n+m,1);
    [~,kvecyuan(i),tvecyuan(i),evecyuan(i)] = Yuan_Varios_Fun(x0,n,m,A1,b1,c1,A2,b2,c2,mystr);
    
    asalvarstr = strcat('.\Salvar\resultsquadratica_yuan',num2str(n),num2str(m),nustr);
    save(asalvarstr,'kvec','tvec','evec','kvecjacobi','tvecjacobi','evecjacobi','kvecintern','tvecintern','evecintern','kvecdois','tvecdois','evecdois','kvecyuan','tvecyuan','evecyuan');
end


clear


%%%%%%%%%% 10x10


%%%%%%%%%quartico %%%%%%%%

alfa = 1/2; 
tol = 1e-4; kmax=1000;
minvaltol = 1;
maxvaltol=1e10;
%30 problemas de n=m=10
n=10; m=10;
mystr = 'quartico';
x0=ones(n+m,1);
kvec=zeros(30,1); %guardar iter
tvec=zeros(30,1); %guardar t
evec=zeros(30,1); %guardar erro

kvecjacobi=zeros(30,1); kvecintern=kvec; kvecdois=kvec; kvecyuan = kvec;
tvecjacobi=zeros(30,1); tvecintern=tvec; tvecdois=tvec; tvecyuan = tvec;
evecjacobi=zeros(30,1); evecintern=evec; evecdois=evec; evecyuan = evec; 

A1 = [];
A2 = [];
                
for i=1:30
    nustr = num2str(i);   
    b1 = (-1.5+i/10)*ones(n,1); b2=b1;
    c1 = (1.5-i/10)*ones(m,1); c2=c1;
    flaghess = 'ifdois'; multhess = 1;
    Geral
    kvec(i) = k;
    tvec(i) = maint;
    evec(i) = erro;
    
    flaghess = 'jacobiintern'; x0=ones(n+m,1);
    Geral
    kvecintern(i) = k;
    tvecintern(i) = maint;
    evecintern(i) = erro;
    
    flaghess = 'jacobidois'; x0=ones(n+m,1);
    Geral
    kvecdois(i) = k;
    tvecdois(i) = maint;
    evecdois(i) = erro;
    
    x0=ones(n+m,1);
    [~,kvecyuan(i),tvecyuan(i),evecyuan(i)] = Yuan_Varios_Fun(x0,n,m,A1,b1,c1,A2,b2,c2,mystr);
    
    asalvarstr = strcat('.\Salvar\resultsquartico',num2str(n),num2str(m),nustr);
    save(asalvarstr,'kvec','tvec','evec','kvecjacobi','tvecjacobi','evecjacobi','kvecintern','tvecintern','evecintern','kvecdois','tvecdois','evecdois','kvecyuan','tvecyuan','evecyuan');
end


clear


%%%%%%% testepaper%%%%%%

alfa = 1/2; 
tol = 1e-4; kmax=1000;
minvaltol = 1;
maxvaltol=1e10;
%30 problemas de n=m=10
n=10; m=10;
mystr = 'testepaper';
x0=ones(n+m,1);
kvec=zeros(30,1); %guardar iter
tvec=zeros(30,1); %guardar t
evec=zeros(30,1); %guardar erro

kvecjacobi=zeros(30,1); kvecintern=kvec; kvecdois=kvec; kvecyuan = kvec;
tvecjacobi=zeros(30,1); tvecintern=tvec; tvecdois=tvec; tvecyuan = tvec;
evecjacobi=zeros(30,1); evecintern=evec; evecdois=evec; evecyuan = evec;

c1=zeros(n,1);
c2=zeros(m,1);
for i=1:30
    nustr = num2str(i);
    fistr = strcat('.\Salvar\matrizestestepaper',num2str(n),num2str(m),nustr);
    quadfile = matfile(fistr);
    A1 = quadfile.A1;
    A2 = quadfile.A2;
    b1 = quadfile.b1;
    b2 = quadfile.b2;
    
    flaghess = 'ifdois'; multhess = 1;
    Geral
    kvec(i) = k;
    tvec(i) = maint;
    evec(i) = erro;
    
    kvecjacobi(i) = kjacobi;
    tvecjacobi(i) = tjacobi;
    evecjacobi(i) = ejacobi;
    
    flaghess = 'jacobiintern'; x0=ones(n+m,1);
    Geral
    kvecintern(i) = k;
    tvecintern(i) = maint;
    evecintern(i) = erro;
    
    flaghess = 'jacobidois'; x0=ones(n+m,1);
    Geral
    kvecdois(i) = k;
    tvecdois(i) = maint;
    evecdois(i) = erro;
    
    x0=ones(n+m,1);
    [~,kvecyuan(i),tvecyuan(i),evecyuan(i)] = Yuan_Varios_Fun(x0,n,m,A1,b1,c1,A2,b2,c2,mystr);
    
    asalvarstr = strcat('.\Salvar\resultstestepaper',num2str(n),num2str(m),nustr);
    save(asalvarstr,'kvec','tvec','evec','kvecjacobi','tvecjacobi','evecjacobi','kvecintern','tvecintern','evecintern','kvecdois','tvecdois','evecdois','kvecyuan','tvecyuan','evecyuan');
end


clear


erroaqui

%%%%%%% vacina %%%%%

alfa = 1/2; 
tol = 1e-4; kmax=1000;
minvaltol = 1;
maxvaltol=1e10;
%30 problemas de n=m=10
n=10; m=10;
mystr = 'vacina';
x0=ones(n+m,1);
kvec=zeros(30,1); %guardar iter
tvec=zeros(30,1); %guardar t
evec=zeros(30,1); %guardar erro

kvecjacobi=zeros(30,1); kvecintern=kvec; kvecdois=kvec; kvecyuan = kvec;
tvecjacobi=zeros(30,1); tvecintern=tvec; tvecdois=tvec; tvecyuan = tvec;
evecjacobi=zeros(30,1); evecintern=evec; evecdois=evec; evecyuan = evec;

A1 = [];
A2 = [];
b1 = 0.45*eye(n,m);
b2 = 0.45*eye(m,n);
                
for i=1:30
    nustr = num2str(i);
    c1 = (1.5-i/10)*ones(n,1);
    c2 = (1.4-i/10)*ones(m,1);
    flaghess = 'ifdois'; multhess = 1;
    Geral
    kvec(i) = k;
    tvec(i) = maint;
    evec(i) = erro;
    
    kvecjacobi(i) = kjacobi;
    tvecjacobi(i) = tjacobi;
    evecjacobi(i) = ejacobi;
    
    flaghess = 'jacobiintern'; x0=ones(n+m,1);
    Geral
    kvecintern(i) = k;
    tvecintern(i) = maint;
    evecintern(i) = erro;
    
    flaghess = 'jacobidois'; x0=ones(n+m,1);
    Geral
    kvecdois(i) = k;
    tvecdois(i) = maint;
    evecdois(i) = erro;
    
    x0=ones(n+m,1);
    [~,kvecyuan(i),tvecyuan(i),evecyuan(i)] = Yuan_Varios_Fun(x0,n,m,A1,b1,c1,A2,b2,c2,mystr);
    
    asalvarstr = strcat('.\Salvar\resultsvacina',num2str(n),num2str(m),nustr);
    save(asalvarstr,'kvec','tvec','evec','kvecjacobi','tvecjacobi','evecjacobi','kvecintern','tvecintern','evecintern','kvecdois','tvecdois','evecdois','kvecyuan','tvecyuan','evecyuan');
end


clear


%%%%%%%%%%%%%%%%%%%%% 50x50 %%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%problemayuan%%%%%%

alfa = 1/2; 
tol = 1e-4; kmax=1000;
minvaltol = 1;
maxvaltol=1e10;
%30 problemas de n=m=10
n=50; m=50;
mystr = 'problemayuan';
x0=ones(n+m,1);
kvec=zeros(30,1); %guardar iter
tvec=zeros(30,1); %guardar t
evec=zeros(30,1); %guardar erro

kvecjacobi=zeros(30,1); kvecintern=kvec; kvecdois=kvec; kvecyuan = kvec;
tvecjacobi=zeros(30,1); tvecintern=tvec; tvecdois=tvec; tvecyuan = tvec;
evecjacobi=zeros(30,1); evecintern=evec; evecdois=evec; evecyuan = evec; 

A1 = zeros(n,n);
A2 = zeros(m,m);
b1 = eye(n,m);
b2 = -eye(m,n);
                
for i=1:30
    nustr = num2str(i);   
    c1 = (-1.5+i/10)*ones(n,1);
    c2 = (1.4-i/10)*ones(m,1);
    flaghess = 'ifdois'; multhess = 1;
    Geral
    kvec(i) = k;
    tvec(i) = maint;
    evec(i) = erro;
    
    kvecjacobi(i) = kjacobi;
    tvecjacobi(i) = tjacobi;
    evecjacobi(i) = ejacobi;
    
    flaghess = 'jacobiintern'; x0=ones(n+m,1);
    Geral
    kvecintern(i) = k;
    tvecintern(i) = maint;
    evecintern(i) = erro;
    
    flaghess = 'jacobidois'; x0=ones(n+m,1);
    Geral
    kvecdois(i) = k;
    tvecdois(i) = maint;
    evecdois(i) = erro;
    
    x0=ones(n+m,1);
    [~,kvecyuan(i),tvecyuan(i),evecyuan(i)] = Yuan_Varios_Fun(x0,n,m,A1,b1,c1,A2,b2,c2,mystr);
    
    asalvarstr = strcat('.\Salvar\resultsproblemayuan',num2str(n),num2str(m),nustr);
    save(asalvarstr,'kvec','tvec','evec','kvecjacobi','tvecjacobi','evecjacobi','kvecintern','tvecintern','evecintern','kvecdois','tvecdois','evecdois','kvecyuan','tvecyuan','evecyuan');
end


clear


%%%quartico %%%%%%%

alfa = 1/2; 
tol = 1e-4; kmax=1000;
minvaltol = 1;
maxvaltol=1e10;
%30 problemas de n=m=10
n=50; m=50;
mystr = 'quartico';
x0=ones(n+m,1);
kvec=zeros(30,1); %guardar iter
tvec=zeros(30,1); %guardar t
evec=zeros(30,1); %guardar erro

kvecjacobi=zeros(30,1); kvecintern=kvec; kvecdois=kvec; kvecyuan = kvec;
tvecjacobi=zeros(30,1); tvecintern=tvec; tvecdois=tvec; tvecyuan = tvec;
evecjacobi=zeros(30,1); evecintern=evec; evecdois=evec; evecyuan = evec; 

A1 = [];
A2 = [];
                
for i=1:30
    nustr = num2str(i);   
    b1 = (-1.5+i/10)*ones(n,1); b2=b1;
    c1 = (1.5-i/10)*ones(m,1); c2=c1;
    flaghess = 'ifdois'; multhess = 1;
    Geral
    kvec(i) = k;
    tvec(i) = maint;
    evec(i) = erro;
    
    
    flaghess = 'jacobiintern'; x0=ones(n+m,1);
    Geral
    kvecintern(i) = k;
    tvecintern(i) = maint;
    evecintern(i) = erro;
    
    flaghess = 'jacobidois'; x0=ones(n+m,1);
    Geral
    kvecdois(i) = k;
    tvecdois(i) = maint;
    evecdois(i) = erro;
    
    x0=ones(n+m,1);
    [~,kvecyuan(i),tvecyuan(i),evecyuan(i)] = Yuan_Varios_Fun(x0,n,m,A1,b1,c1,A2,b2,c2,mystr);
    
    asalvarstr = strcat('.\Salvar\resultsquartico',num2str(n),num2str(m),nustr);
    save(asalvarstr,'kvec','tvec','evec','kvecjacobi','tvecjacobi','evecjacobi','kvecintern','tvecintern','evecintern','kvecdois','tvecdois','evecdois','kvecyuan','tvecyuan','evecyuan');
end


clear


%%%%%%%% quadratica yuan %%%%%%%


alfa = 1/2; 
tol = 1e-4; kmax=1000;
minvaltol = 1;
maxvaltol=1e10;
%30 problemas de n=m=10
n=50; m=50;
mystr = 'quadratica_yuan';
x0=ones(n+m,1);
kvec=zeros(30,1); %guardar iter
tvec=zeros(30,1); %guardar t
evec=zeros(30,1); %guardar erro

kvecjacobi=zeros(30,1); kvecintern=kvec; kvecdois=kvec; kvecyuan = kvec;
tvecjacobi=zeros(30,1); tvecintern=tvec; tvecdois=tvec; tvecyuan = tvec;
evecjacobi=zeros(30,1); evecintern=evec; evecdois=evec; evecyuan = evec;

A1 = [];
A2 = [];
                
for i=1:30
    nustr = num2str(i);   
    b1 = (-1.5+i/10)*ones(n,1); b2=b1;
    c1 = (1.5-i/10)*ones(m,1); c2=c1;
    flaghess = 'ifdois'; multhess = 1;
    Geral
    kvec(i) = k;
    tvec(i) = maint;
    evec(i) = erro;
    
    kvecjacobi(i) = kjacobi;
    tvecjacobi(i) = tjacobi;
    evecjacobi(i) = ejacobi;
    
    flaghess = 'jacobiintern'; x0=ones(n+m,1);
    Geral
    kvecintern(i) = k;
    tvecintern(i) = maint;
    evecintern(i) = erro;
    
    flaghess = 'jacobidois'; x0=ones(n+m,1);
    Geral
    kvecdois(i) = k;
    tvecdois(i) = maint;
    evecdois(i) = erro;
    
    x0=ones(n+m,1);
    [~,kvecyuan(i),tvecyuan(i),evecyuan(i)] = Yuan_Varios_Fun(x0,n,m,A1,b1,c1,A2,b2,c2,mystr);
    
    asalvarstr = strcat('.\Salvar\resultsquadratica_yuan',num2str(n),num2str(m),nustr);
    save(asalvarstr,'kvec','tvec','evec','kvecjacobi','tvecjacobi','evecjacobi','kvecintern','tvecintern','evecintern','kvecdois','tvecdois','evecdois','kvecyuan','tvecyuan','evecyuan');
end


clear

%%%%%%%% testepaper%%%%


alfa = 1/2; 
tol = 1e-4; kmax=1000;
minvaltol = 1;
maxvaltol=1e10;
%30 problemas de n=m=10
n=50; m=50;
mystr = 'testepaper';
x0=ones(n+m,1);
kvec=zeros(30,1); %guardar iter
tvec=zeros(30,1); %guardar t
evec=zeros(30,1); %guardar erro

kvecjacobi=zeros(30,1); kvecintern=kvec; kvecdois=kvec; kvecyuan = kvec;
tvecjacobi=zeros(30,1); tvecintern=tvec; tvecdois=tvec; tvecyuan = tvec;
evecjacobi=zeros(30,1); evecintern=evec; evecdois=evec; evecyuan = evec; 

c1=zeros(n,1);
c2=zeros(m,1);
for i=1:30
    nustr = num2str(i);
    fistr = strcat('.\Salvar\matrizestestepaper',num2str(n),num2str(m),nustr);
    quadfile = matfile(fistr);
    A1 = quadfile.A1;
    A2 = quadfile.A2;
    b1 = quadfile.b1;
    b2 = quadfile.b2;
    
    flaghess = 'ifdois'; multhess = 1;
    Geral
    kvec(i) = k;
    tvec(i) = maint;
    evec(i) = erro;
    
    kvecjacobi(i) = kjacobi;
    tvecjacobi(i) = tjacobi;
    evecjacobi(i) = ejacobi;
    
    flaghess = 'jacobiintern'; x0=ones(n+m,1);
    Geral
    kvecintern(i) = k;
    tvecintern(i) = maint;
    evecintern(i) = erro;
    
    flaghess = 'jacobidois'; x0=ones(n+m,1);
    Geral
    kvecdois(i) = k;
    tvecdois(i) = maint;
    evecdois(i) = erro;
    
    x0=ones(n+m,1);
    [~,kvecyuan(i),tvecyuan(i),evecyuan(i)] = Yuan_Varios_Fun(x0,n,m,A1,b1,c1,A2,b2,c2,mystr);
    
    asalvarstr = strcat('.\Salvar\resultstestepaper',num2str(n),num2str(m),nustr);
    save(asalvarstr,'kvec','tvec','evec','kvecjacobi','tvecjacobi','evecjacobi','kvecintern','tvecintern','evecintern','kvecdois','tvecdois','evecdois','kvecyuan','tvecyuan','evecyuan');
end    



clear

%%%%%%%% vacina %%%%%%%%%


alfa = 1/2; 
tol = 1e-4; kmax=1000;
minvaltol = 1;
maxvaltol=1e10;
%30 problemas de n=m=10
n=50; m=50;
mystr = 'vacina';
x0=ones(n+m,1);
kvec=zeros(30,1); %guardar iter
tvec=zeros(30,1); %guardar t
evec=zeros(30,1); %guardar erro

kvecjacobi=zeros(30,1); kvecintern=kvec; kvecdois=kvec; kvecyuan = kvec;
tvecjacobi=zeros(30,1); tvecintern=tvec; tvecdois=tvec; tvecyuan = tvec;
evecjacobi=zeros(30,1); evecintern=evec; evecdois=evec; evecyuan = evec;

A1 = [];
A2 = [];
b1 = 0.45*eye(n,m);
b2 = 0.45*eye(m,n);
                
for i=1:30
    nustr = num2str(i);
    c1 = (1.5-i/10)*ones(n,1);
    c2 = (1.4-i/10)*ones(m,1);
    flaghess = 'ifdois'; multhess = 1;
    Geral
    kvec(i) = k;
    tvec(i) = maint;
    evec(i) = erro;
    
    kvecjacobi(i) = kjacobi;
    tvecjacobi(i) = tjacobi;
    evecjacobi(i) = ejacobi;
    
    flaghess = 'jacobiintern'; x0=ones(n+m,1);
    Geral
    kvecintern(i) = k;
    tvecintern(i) = maint;
    evecintern(i) = erro;
    
    flaghess = 'jacobidois'; x0=ones(n+m,1);
    Geral
    kvecdois(i) = k;
    tvecdois(i) = maint;
    evecdois(i) = erro;
    
    x0=ones(n+m,1);
    [~,kvecyuan(i),tvecyuan(i),evecyuan(i)] = Yuan_Varios_Fun(x0,n,m,A1,b1,c1,A2,b2,c2,mystr);
    
    asalvarstr = strcat('.\Salvar\resultsvacina',num2str(n),num2str(m),nustr);
    save(asalvarstr,'kvec','tvec','evec','kvecjacobi','tvecjacobi','evecjacobi','kvecintern','tvecintern','evecintern','kvecdois','tvecdois','evecdois','kvecyuan','tvecyuan','evecyuan');
end


clear


%%%%%%%%%%%%%%%%%% 100x100 %%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%% problema yuan %%%%%%%

alfa = 1/2; 
tol = 1e-4; kmax=1000;
minvaltol = 1;
maxvaltol=1e10;
%30 problemas de n=m=10
n=100; m=100;
mystr = 'problemayuan';
x0=ones(n+m,1);
kvec=zeros(30,1); %guardar iter
tvec=zeros(30,1); %guardar t
evec=zeros(30,1); %guardar erro

kvecjacobi=zeros(30,1); kvecintern=kvec; kvecdois=kvec; kvecyuan = kvec;
tvecjacobi=zeros(30,1); tvecintern=tvec; tvecdois=tvec; tvecyuan = tvec;
evecjacobi=zeros(30,1); evecintern=evec; evecdois=evec; evecyuan = evec; 

A1 = zeros(n,n);
A2 = zeros(m,m);
b1 = eye(n,m);
b2 = -eye(m,n);
                
for i=1:30
    nustr = num2str(i);   
    c1 = (-1.5+i/10)*ones(n,1);
    c2 = (1.4-i/10)*ones(m,1);
    flaghess = 'ifdois'; multhess = 1;
    Geral
    kvec(i) = k;
    tvec(i) = maint;
    evec(i) = erro;
    
    kvecjacobi(i) = kjacobi;
    tvecjacobi(i) = tjacobi;
    evecjacobi(i) = ejacobi;
    
    flaghess = 'jacobiintern'; x0=ones(n+m,1);
    Geral
    kvecintern(i) = k;
    tvecintern(i) = maint;
    evecintern(i) = erro;
    
    flaghess = 'jacobidois'; x0=ones(n+m,1);
    Geral
    kvecdois(i) = k;
    tvecdois(i) = maint;
    evecdois(i) = erro;
    
    x0=ones(n+m,1);
    [~,kvecyuan(i),tvecyuan(i),evecyuan(i)] = Yuan_Varios_Fun(x0,n,m,A1,b1,c1,A2,b2,c2,mystr);
    
    asalvarstr = strcat('.\Salvar\resultsproblemayuan',num2str(n),num2str(m),nustr);
    save(asalvarstr,'kvec','tvec','evec','kvecjacobi','tvecjacobi','evecjacobi','kvecintern','tvecintern','evecintern','kvecdois','tvecdois','evecdois','kvecyuan','tvecyuan','evecyuan');
end


clear

%%%%%%% quartico %%%%%%

alfa = 1/2; 
tol = 1e-4; kmax=1000;
minvaltol = 1;
maxvaltol=1e10;
%30 problemas de n=m=10
n=100; m=100;
mystr = 'quartico';
x0=ones(n+m,1);
kvec=zeros(30,1); %guardar iter
tvec=zeros(30,1); %guardar t
evec=zeros(30,1); %guardar erro

kvecjacobi=zeros(30,1); kvecintern=kvec; kvecdois=kvec; kvecyuan = kvec;
tvecjacobi=zeros(30,1); tvecintern=tvec; tvecdois=tvec; tvecyuan = tvec;
evecjacobi=zeros(30,1); evecintern=evec; evecdois=evec; evecyuan = evec; 

A1 = [];
A2 = [];
                
for i=1:30
    nustr = num2str(i);   
    b1 = (-1.5+i/10)*ones(n,1); b2=b1;
    c1 = (1.5-i/10)*ones(m,1); c2=c1;
    flaghess = 'ifdois'; multhess = 1;
    Geral
    kvec(i) = k;
    tvec(i) = maint;
    evec(i) = erro;
    
    flaghess = 'jacobiintern'; x0=ones(n+m,1);
    Geral
    kvecintern(i) = k;
    tvecintern(i) = maint;
    evecintern(i) = erro;
    
    flaghess = 'jacobidois'; x0=ones(n+m,1);
    Geral
    kvecdois(i) = k;
    tvecdois(i) = maint;
    evecdois(i) = erro;
    
    x0=ones(n+m,1);
    [~,kvecyuan(i),tvecyuan(i),evecyuan(i)] = Yuan_Varios_Fun(x0,n,m,A1,b1,c1,A2,b2,c2,mystr);
    
    asalvarstr = strcat('.\Salvar\resultsquartico',num2str(n),num2str(m),nustr);
    save(asalvarstr,'kvec','tvec','evec','kvecjacobi','tvecjacobi','evecjacobi','kvecintern','tvecintern','evecintern','kvecdois','tvecdois','evecdois','kvecyuan','tvecyuan','evecyuan');
end



clear

%%%%% quadratica yuan %%%%%%


alfa = 1/2; 
tol = 1e-4; kmax=1000;
minvaltol = 1;
maxvaltol=1e10;
%30 problemas de n=m=10
n=100; m=100;
mystr = 'quadratica_yuan';
x0=ones(n+m,1);
kvec=zeros(30,1); %guardar iter
tvec=zeros(30,1); %guardar t
evec=zeros(30,1); %guardar erro

kvecjacobi=zeros(30,1); kvecintern=kvec; kvecdois=kvec; kvecyuan = kvec;
tvecjacobi=zeros(30,1); tvecintern=tvec; tvecdois=tvec; tvecyuan = tvec;
evecjacobi=zeros(30,1); evecintern=evec; evecdois=evec; evecyuan = evec; 

A1 = [];
A2 = [];
                
for i=1:30
    nustr = num2str(i);   
    b1 = (-1.5+i/10)*ones(n,1); b2=b1;
    c1 = (1.5-i/10)*ones(m,1); c2=c1;
    flaghess = 'ifdois'; multhess = 1;
    Geral
    kvec(i) = k;
    tvec(i) = maint;
    evec(i) = erro;
    
    kvecjacobi(i) = kjacobi;
    tvecjacobi(i) = tjacobi;
    evecjacobi(i) = ejacobi;
    
    flaghess = 'jacobiintern'; x0=ones(n+m,1);
    Geral
    kvecintern(i) = k;
    tvecintern(i) = maint;
    evecintern(i) = erro;
    
    flaghess = 'jacobidois'; x0=ones(n+m,1);
    Geral
    kvecdois(i) = k;
    tvecdois(i) = maint;
    evecdois(i) = erro;
    
    x0=ones(n+m,1);
    [~,kvecyuan(i),tvecyuan(i),evecyuan(i)] = Yuan_Varios_Fun(x0,n,m,A1,b1,c1,A2,b2,c2,mystr);
    
    asalvarstr = strcat('.\Salvar\resultsquadratica_yuan',num2str(n),num2str(m),nustr);
    save(asalvarstr,'kvec','tvec','evec','kvecjacobi','tvecjacobi','evecjacobi','kvecintern','tvecintern','evecintern','kvecdois','tvecdois','evecdois','kvecyuan','tvecyuan','evecyuan');
end


clear

%%%%%% testepaper%%%%%


alfa = 1/2; 
tol = 1e-4; kmax=1000;
minvaltol = 1;
maxvaltol=1e10;
%30 problemas de n=m=10
n=100; m=100;
mystr = 'testepaper';
x0=ones(n+m,1);
kvec=zeros(30,1); %guardar iter
tvec=zeros(30,1); %guardar t
evec=zeros(30,1); %guardar erro

kvecjacobi=zeros(30,1); kvecintern=kvec; kvecdois=kvec; kvecyuan = kvec;
tvecjacobi=zeros(30,1); tvecintern=tvec; tvecdois=tvec; tvecyuan = tvec;
evecjacobi=zeros(30,1); evecintern=evec; evecdois=evec; evecyuan = evec;

c1=zeros(n,1);
c2=zeros(m,1);
for i=1:30
    nustr = num2str(i);
    fistr = strcat('.\Salvar\matrizestestepaper',num2str(n),num2str(m),nustr);
    quadfile = matfile(fistr);
    A1 = quadfile.A1;
    A2 = quadfile.A2;
    b1 = quadfile.b1;
    b2 = quadfile.b2;
    
    flaghess = 'ifdois'; multhess = 1;
    Geral
    kvec(i) = k;
    tvec(i) = maint;
    evec(i) = erro;
    
    kvecjacobi(i) = kjacobi;
    tvecjacobi(i) = tjacobi;
    evecjacobi(i) = ejacobi;
    
    flaghess = 'jacobiintern'; x0=ones(n+m,1);
    Geral
    kvecintern(i) = k;
    tvecintern(i) = maint;
    evecintern(i) = erro;
    
    flaghess = 'jacobidois'; x0=ones(n+m,1);
    Geral
    kvecdois(i) = k;
    tvecdois(i) = maint;
    evecdois(i) = erro;
    
    x0=ones(n+m,1);
    [~,kvecyuan(i),tvecyuan(i),evecyuan(i)] = Yuan_Varios_Fun(x0,n,m,A1,b1,c1,A2,b2,c2,mystr);
    
    asalvarstr = strcat('.\Salvar\resultstestepaper',num2str(n),num2str(m),nustr);
    save(asalvarstr,'kvec','tvec','evec','kvecjacobi','tvecjacobi','evecjacobi','kvecintern','tvecintern','evecintern','kvecdois','tvecdois','evecdois','kvecyuan','tvecyuan','evecyuan');
end


clear


%%%%% vacina %%%%%%%%

alfa = 1/2; 
tol = 1e-4; kmax=1000;
minvaltol = 1;
maxvaltol=1e10;
%30 problemas de n=m=10
n=100; m=100;
mystr = 'vacina';
x0=ones(n+m,1);
kvec=zeros(30,1); %guardar iter
tvec=zeros(30,1); %guardar t
evec=zeros(30,1); %guardar erro

kvecjacobi=zeros(30,1); kvecintern=kvec; kvecdois=kvec; kvecyuan = kvec;
tvecjacobi=zeros(30,1); tvecintern=tvec; tvecdois=tvec; tvecyuan = tvec;
evecjacobi=zeros(30,1); evecintern=evec; evecdois=evec; evecyuan = evec;

A1 = [];
A2 = [];
b1 = 0.45*eye(n,m);
b2 = 0.45*eye(m,n);
                
for i=1:30
    nustr = num2str(i);
    c1 = (1.5-i/10)*ones(n,1);
    c2 = (1.4-i/10)*ones(m,1);
    flaghess = 'ifdois'; multhess = 1;
    Geral
    kvec(i) = k;
    tvec(i) = maint;
    evec(i) = erro;
    
    kvecjacobi(i) = kjacobi;
    tvecjacobi(i) = tjacobi;
    evecjacobi(i) = ejacobi;
    
    flaghess = 'jacobiintern'; x0=ones(n+m,1);
    Geral
    kvecintern(i) = k;
    tvecintern(i) = maint;
    evecintern(i) = erro;
    
    flaghess = 'jacobidois'; x0=ones(n+m,1);
    Geral
    kvecdois(i) = k;
    tvecdois(i) = maint;
    evecdois(i) = erro;
    
    x0=ones(n+m,1);
    [~,kvecyuan(i),tvecyuan(i),evecyuan(i)] = Yuan_Varios_Fun(x0,n,m,A1,b1,c1,A2,b2,c2,mystr);
    
    asalvarstr = strcat('.\Salvar\resultsvacina',num2str(n),num2str(m),nustr);
    save(asalvarstr,'kvec','tvec','evec','kvecjacobi','tvecjacobi','evecjacobi','kvecintern','tvecintern','evecintern','kvecdois','tvecdois','evecdois','kvecyuan','tvecyuan','evecyuan');
end


clear



%%%%%%%%%%%%%%%%%% 1x1 %%%%%%%%%%%%%%



%%%%%%% problema yuan %%%%%%%

alfa = 1/2; 
tol = 1e-4; kmax=1000;
minvaltol = 1;
maxvaltol=1e10;
%30 problemas de n=m=10
n=1; m=1;
mystr = 'problemayuan';
x0=ones(n+m,1);
kvec=zeros(30,1); %guardar iter
tvec=zeros(30,1); %guardar t
evec=zeros(30,1); %guardar erro
xvec=zeros(30,2); %guardar ponto (soh qd n=m=1)

kvecjacobi=zeros(30,1); kvecintern=kvec; kvecdois=kvec; kvecyuan = kvec;
tvecjacobi=zeros(30,1); tvecintern=tvec; tvecdois=tvec; tvecyuan = tvec;
evecjacobi=zeros(30,1); evecintern=evec; evecdois=evec; evecyuan = evec;
xvecjacobi=zeros(30,2); xvecintern=xvec; xvecdois=xvec; xvecyuan = xvec;

A1 = zeros(n,n);
A2 = zeros(m,m);
b1 = eye(n,m);
b2 = -eye(m,n);
                
for i=1:30
    nustr = num2str(i);
    c1 = (-1.5+i/10)*ones(n,1);
    c2 = (1.4-i/10)*ones(m,1);
    flaghess = 'ifdois'; multhess = 1;
    Geral
    kvec(i) = k;
    tvec(i) = maint;
    evec(i) = erro;
    xvec(i,1) = x0(1);
    xvec(i,2) = x0(2);
    
    kvecjacobi(i) = kjacobi;
    tvecjacobi(i) = tjacobi;
    evecjacobi(i) = ejacobi;
    xvecjacobi(i,1) = xjacobi(1);
    xvecjacobi(i,2) = xjacobi(2);
    
    flaghess = 'jacobiintern'; x0=ones(n+m,1);
    Geral
    kvecintern(i) = k;
    tvecintern(i) = maint;
    evecintern(i) = erro;
    xvecintern(i,1) = x0(1);
    xvecintern(i,2) = x0(2);
    
    flaghess = 'jacobidois'; x0=ones(n+m,1);
    Geral
    kvecdois(i) = k;
    tvecdois(i) = maint;
    evecdois(i) = erro;
    xvecdois(i,1) = x0(1);
    xvecdois(i,2) = x0(2);
    
    x0=ones(n+m,1);
    [xvecret,kvecyuan(i),tvecyuan(i),evecyuan(i)] = Yuan_Varios_Fun(x0,n,m,A1,b1,c1,A2,b2,c2,mystr);
    xvecyuan(i,1)=xvecret(1);
    xvecyuan(i,2)=xvecret(2);
      
    asalvarstr = strcat('.\Salvar\resultsproblemayuan',num2str(n),num2str(m),nustr);
    save(asalvarstr,'kvec','tvec','evec','xvec','kvecjacobi','tvecjacobi','evecjacobi','xvecjacobi','kvecintern','tvecintern','evecintern','xvecintern','kvecdois','tvecdois','evecdois','xvecdois','kvecyuan','tvecyuan','evecyuan','xvecyuan');
end


clear


%%%%%testepaper  %%%

alfa = 1/2; 
tol = 1e-4; kmax=1000;
minvaltol = 1;
maxvaltol=1e10;
%primeiro 30 problemas de n=m=1
n=1; m=1;
x0 = ones(n+m,1);
mystr='testepaper';
kvec=zeros(30,1); %guardar iter
tvec=zeros(30,1); %guardar t
evec=zeros(30,1); %guardar erro
xvec=zeros(30,2); %guardar ponto (soh qd n=m=1)

kvecjacobi=zeros(30,1); kvecintern=kvec; kvecdois=kvec; kvecyuan = kvec;
tvecjacobi=zeros(30,1); tvecintern=tvec; tvecdois=tvec; tvecyuan = tvec;
evecjacobi=zeros(30,1); evecintern=evec; evecdois=evec; evecyuan = evec;
xvecjacobi=zeros(30,2); xvecintern=xvec; xvecdois=xvec; xvecyuan = xvec;

c1=zeros(n,1);
c2=zeros(m,1);
for i=1:30
    nustr = num2str(i);
    fistr = strcat('.\Salvar\matrizestestepaper',num2str(n),num2str(m),nustr);
    quadfile = matfile(fistr);
    A1 = quadfile.A1;
    A2 = quadfile.A2;
    b1 = quadfile.b1;
    b2 = quadfile.b2;
    
    flaghess = 'ifdois'; multhess = 1;
    Geral
    kvec(i) = k;
    tvec(i) = maint;
    evec(i) = erro;
    xvec(i,1) = x0(1);
    xvec(i,2) = x0(2);
    
    kvecjacobi(i) = kjacobi;
    tvecjacobi(i) = tjacobi;
    evecjacobi(i) = ejacobi;
    xvecjacobi(i,1) = xjacobi(1);
    xvecjacobi(i,2) = xjacobi(2);
    
    flaghess = 'jacobiintern'; x0=ones(n+m,1);
    Geral
    kvecintern(i) = k;
    tvecintern(i) = maint;
    evecintern(i) = erro;
    xvecintern(i,1) = x0(1);
    xvecintern(i,2) = x0(2);
    
    flaghess = 'jacobidois'; x0=ones(n+m,1);
    Geral
    kvecdois(i) = k;
    tvecdois(i) = maint;
    evecdois(i) = erro;
    xvecdois(i,1) = x0(1);
    xvecdois(i,2) = x0(2);
    
    x0=ones(n+m,1);
    [xvecret,kvecyuan(i),tvecyuan(i),evecyuan(i)] = Yuan_Varios_Fun(x0,n,m,A1,b1,c1,A2,b2,c2,mystr);
    xvecyuan(i,1)=xvecret(1);
    xvecyuan(i,2)=xvecret(2);
      
    asalvarstr = strcat('.\Salvar\resultstestepaper',num2str(n),num2str(m),nustr);
    save(asalvarstr,'kvec','tvec','evec','xvec','kvecjacobi','tvecjacobi','evecjacobi','xvecjacobi','kvecintern','tvecintern','evecintern','xvecintern','kvecdois','tvecdois','evecdois','xvecdois','kvecyuan','tvecyuan','evecyuan','xvecyuan');
end


clear

%%%%%%%% Vacina %%%%%%


alfa = 1/2; 
tol = 1e-4; kmax=1000;
minvaltol = 1;
maxvaltol=1e10;
%30 problemas de n=m=10
n=1; m=1;
mystr = 'vacina';
x0=ones(n+m,1);
kvec=zeros(30,1); %guardar iter
tvec=zeros(30,1); %guardar t
evec=zeros(30,1); %guardar erro
xvec=zeros(30,2); %guardar ponto (soh qd n=m=1)

kvecjacobi=zeros(30,1); kvecintern=kvec; kvecdois=kvec; kvecyuan = kvec;
tvecjacobi=zeros(30,1); tvecintern=tvec; tvecdois=tvec; tvecyuan = tvec;
evecjacobi=zeros(30,1); evecintern=evec; evecdois=evec; evecyuan = evec;
xvecjacobi=zeros(30,2); xvecintern=xvec; xvecdois=xvec; xvecyuan = xvec;

A1 = [];
A2 = [];
b1 = 0.45*eye(n,m);
b2 = 0.45*eye(m,n);
                
for i=1:30
    nustr = num2str(i);
    c1 = (1.5-i/10)*ones(n,1);
    c2 = (1.4-i/10)*ones(m,1);
    flaghess = 'ifdois'; multhess = 1;
    Geral
    kvec(i) = k;
    tvec(i) = maint;
    evec(i) = erro;
    xvec(i,1) = x0(1);
    xvec(i,2) = x0(2);
    
    kvecjacobi(i) = kjacobi;
    tvecjacobi(i) = tjacobi;
    evecjacobi(i) = ejacobi;
    xvecjacobi(i,1) = xjacobi(1);
    xvecjacobi(i,2) = xjacobi(2);
    
    flaghess = 'jacobiintern'; x0=ones(n+m,1);
    Geral
    kvecintern(i) = k;
    tvecintern(i) = maint;
    evecintern(i) = erro;
    xvecintern(i,1) = x0(1);
    xvecintern(i,2) = x0(2);
    
    flaghess = 'jacobidois'; x0=ones(n+m,1);
    Geral
    kvecdois(i) = k;
    tvecdois(i) = maint;
    evecdois(i) = erro;
    xvecdois(i,1) = x0(1);
    xvecdois(i,2) = x0(2);
    
    x0=ones(n+m,1);
    [xvecret,kvecyuan(i),tvecyuan(i),evecyuan(i)] = Yuan_Varios_Fun(x0,n,m,A1,b1,c1,A2,b2,c2,mystr);
    xvecyuan(i,1)=xvecret(1);
    xvecyuan(i,2)=xvecret(2);
      
    asalvarstr = strcat('.\Salvar\resultsvacina',num2str(n),num2str(m),nustr);
    save(asalvarstr,'kvec','tvec','evec','xvec','kvecjacobi','tvecjacobi','evecjacobi','xvecjacobi','kvecintern','tvecintern','evecintern','xvecintern','kvecdois','tvecdois','evecdois','xvecdois','kvecyuan','tvecyuan','evecyuan','xvecyuan');
end


clear

%%%%%%%%%%%% 10x10 %%%%%%%%%%%%%%%%%%%%


alfa = 1/2; 
tol = 1e-4; kmax=1000;
minvaltol = 1;
maxvaltol=1e10;
%30 problemas de n=m=10
n=10; m=10;
mystr = 'problemayuan';
x0=ones(n+m,1);
kvec=zeros(30,1); %guardar iter
tvec=zeros(30,1); %guardar t
evec=zeros(30,1); %guardar erro

kvecjacobi=zeros(30,1); kvecintern=kvec; kvecdois=kvec; kvecyuan = kvec;
tvecjacobi=zeros(30,1); tvecintern=tvec; tvecdois=tvec; tvecyuan = tvec;
evecjacobi=zeros(30,1); evecintern=evec; evecdois=evec; evecyuan = evec; 

A1 = zeros(n,n);
A2 = zeros(m,m);
b1 = eye(n,m);
b2 = -eye(m,n);
                
for i=1:30
    nustr = num2str(i);   
    c1 = (-1.5+i/10)*ones(n,1);
    c2 = (1.4-i/10)*ones(m,1);
    flaghess = 'ifdois'; multhess = 1;
    Geral
    kvec(i) = k;
    tvec(i) = maint;
    evec(i) = erro;
    
    kvecjacobi(i) = kjacobi;
    tvecjacobi(i) = tjacobi;
    evecjacobi(i) = ejacobi;
    
    flaghess = 'jacobiintern'; x0=ones(n+m,1);
    Geral
    kvecintern(i) = k;
    tvecintern(i) = maint;
    evecintern(i) = erro;
    
    flaghess = 'jacobidois'; x0=ones(n+m,1);
    Geral
    kvecdois(i) = k;
    tvecdois(i) = maint;
    evecdois(i) = erro;
    
    x0=ones(n+m,1);
    [~,kvecyuan(i),tvecyuan(i),evecyuan(i)] = Yuan_Varios_Fun(x0,n,m,A1,b1,c1,A2,b2,c2,mystr);
    
    asalvarstr = strcat('.\Salvar\resultsproblemayuan',num2str(n),num2str(m),nustr);
    save(asalvarstr,'kvec','tvec','evec','kvecjacobi','tvecjacobi','evecjacobi','kvecintern','tvecintern','evecintern','kvecdois','tvecdois','evecdois','kvecyuan','tvecyuan','evecyuan');
end


clear



%%%%%%%%%%%
%cubico
%%%%%%%%%%%


alfa = 1/2; 
tol = 1e-4; kmax=1000;
minvaltol = 1;
maxvaltol=1e10;
%30 problemas de n=m=10
n=1; m=1;
mystr = 'cubico';
x0 = ones(n+m,1);

A1 = [];
A2 = [];
b1 = [];
b2 = [];
c1=[];
c2=[];
                
       
    flaghess = 'ifdois'; multhess = 1;
    Geral
    kvec = k;
    tvec = maint;
    evec = erro;
    
    kvecjacobi = kjacobi;
    tvecjacobi = tjacobi;
    evecjacobi = ejacobi;
    
    flaghess = 'jacobiintern'; x0=ones(n+m,1);
    Geral
    kvecintern = k;
    tvecintern = maint;
    evecintern = erro;
    
    flaghess = 'jacobidois'; x0=ones(n+m,1);
    Geral
    kvecdois = k;
    tvecdois = maint;
    evecdois = erro;
    
    x0=ones(n+m,1);
    [~,kvecyuan,tvecyuan,evecyuan] = Yuan_Varios_Fun(x0,n,m,A1,b1,c1,A2,b2,c2,mystr);
    
    asalvarstr = strcat('.\Salvar\resultscubico',num2str(n),num2str(m));
    save(asalvarstr,'kvec','tvec','evec','kvecjacobi','tvecjacobi','evecjacobi','kvecintern','tvecintern','evecintern','kvecdois','tvecdois','evecdois','kvecyuan','tvecyuan','evecyuan');


clear


alfa = 1/2; 
tol = 1e-4; kmax=1000;
minvaltol = 1;
maxvaltol=1e10;
%30 problemas de n=m=10
n=10; m=10;
x0 = ones(n+m,1);
mystr = 'cubico';

A1 = [];
A2 = [];
b1 = [];
b2 = [];
c1=[];
c2=[];
                
       
    flaghess = 'ifdois'; multhess = 1;
    Geral
    kvec = k;
    tvec = maint;
    evec = erro;
    
    kvecjacobi = kjacobi;
    tvecjacobi = tjacobi;
    evecjacobi = ejacobi;
    
    flaghess = 'jacobiintern'; x0=ones(n+m,1);
    Geral
    kvecintern = k;
    tvecintern = maint;
    evecintern = erro;
    
    flaghess = 'jacobidois'; x0=ones(n+m,1);
    Geral
    kvecdois = k;
    tvecdois = maint;
    evecdois = erro;
    
    x0=ones(n+m,1);
    [~,kvecyuan,tvecyuan,evecyuan] = Yuan_Varios_Fun(x0,n,m,A1,b1,c1,A2,b2,c2,mystr);
    
    asalvarstr = strcat('.\Salvar\resultscubico',num2str(n),num2str(m));
    save(asalvarstr,'kvec','tvec','evec','kvecjacobi','tvecjacobi','evecjacobi','kvecintern','tvecintern','evecintern','kvecdois','tvecdois','evecdois','kvecyuan','tvecyuan','evecyuan');


clear


alfa = 1/2; 
tol = 1e-4; kmax=1000;
minvaltol = 1;
maxvaltol=1e10;
%30 problemas de n=m=10
n=50; m=50;
mystr = 'cubico';
x0 = ones(n+m,1);

A1 = [];
A2 = [];
b1 = [];
b2 = [];
c1=[];
c2=[];
                
       
    flaghess = 'ifdois'; multhess = 1;
    Geral
    kvec = k;
    tvec = maint;
    evec = erro;
    
    kvecjacobi = kjacobi;
    tvecjacobi = tjacobi;
    evecjacobi = ejacobi;
    
    flaghess = 'jacobiintern'; x0=ones(n+m,1);
    Geral
    kvecintern = k;
    tvecintern = maint;
    evecintern = erro;
    
    flaghess = 'jacobidois'; x0=ones(n+m,1);
    Geral
    kvecdois = k;
    tvecdois = maint;
    evecdois = erro;
    
    x0=ones(n+m,1);
   [~,kvecyuan,tvecyuan,evecyuan] = Yuan_Varios_Fun(x0,n,m,A1,b1,c1,A2,b2,c2,mystr);
    
    asalvarstr = strcat('.\Salvar\resultscubico',num2str(n),num2str(m));
    save(asalvarstr,'kvec','tvec','evec','kvecjacobi','tvecjacobi','evecjacobi','kvecintern','tvecintern','evecintern','kvecdois','tvecdois','evecdois','kvecyuan','tvecyuan','evecyuan');



clear


alfa = 1/2; 
tol = 1e-4; kmax=1000;
minvaltol = 1;
maxvaltol=1e10;
%30 problemas de n=m=10
n=100; m=100;
mystr = 'cubico';
x0 = ones(n+m,1);

A1 = [];
A2 = [];
b1 = [];
b2 = [];
c1=[];
c2=[];
                
       
    flaghess = 'ifdois'; multhess = 1;
    Geral
    kvec = k;
    tvec = maint;
    evec = erro;
    
    kvecjacobi = kjacobi;
    tvecjacobi = tjacobi;
    evecjacobi = ejacobi;
    
    flaghess = 'jacobiintern'; x0=ones(n+m,1);
    Geral
    kvecintern = k;
    tvecintern = maint;
    evecintern = erro;
    
    flaghess = 'jacobidois'; x0=ones(n+m,1);
    Geral
    kvecdois = k;
    tvecdois = maint;
    evecdois = erro;
    
    x0=ones(n+m,1);
    [~,kvecyuan,tvecyuan,evecyuan] = Yuan_Varios_Fun(x0,n,m,A1,b1,c1,A2,b2,c2,mystr);
    
    asalvarstr = strcat('.\Salvar\resultscubico',num2str(n),num2str(m));
    save(asalvarstr,'kvec','tvec','evec','kvecjacobi','tvecjacobi','evecjacobi','kvecintern','tvecintern','evecintern','kvecdois','tvecdois','evecdois','kvecyuan','tvecyuan','evecyuan');



clear


alfa = 1/2; 
tol = 1e-4; kmax=1000;
minvaltol = 1;
maxvaltol=1e10;
%30 problemas de n=m=10
n=250; m=250;
mystr = 'cubico';
x0 = ones(n+m,1);

A1 = [];
A2 = [];
b1 = [];
b2 = [];
c1=[];
c2=[];
                
       
    flaghess = 'ifdois'; multhess = 1;
    Geral
    kvec = k;
    tvec = maint;
    evec = erro;
    
    kvecjacobi = kjacobi;
    tvecjacobi = tjacobi;
    evecjacobi = ejacobi;
    
    flaghess = 'jacobiintern'; x0=ones(n+m,1);
    Geral
    kvecintern = k;
    tvecintern = maint;
    evecintern = erro;
    
    flaghess = 'jacobidois'; x0=ones(n+m,1);
    Geral
    kvecdois = k;
    tvecdois = maint;
    evecdois = erro;
    
    x0=ones(n+m,1);
    [~,kvecyuan,tvecyuan,evecyuan] = Yuan_Varios_Fun(x0,n,m,A1,b1,c1,A2,b2,c2,mystr);
    
    asalvarstr = strcat('.\Salvar\resultscubico',num2str(n),num2str(m));
    save(asalvarstr,'kvec','tvec','evec','kvecjacobi','tvecjacobi','evecjacobi','kvecintern','tvecintern','evecintern','kvecdois','tvecdois','evecdois','kvecyuan','tvecyuan','evecyuan');



clear



%%%%%%%%%%%%%%%%
%quadtestes1
%%%%%%%%%%%%%%%%


alfa = 1/2; 
tol = 1e-4; kmax=1000;
minvaltol = 1;
maxvaltol=1e10;
%30 problemas de n=m=10
n=1; m=1;
mystr = 'quadtestes1';
x0 = ones(n+m,1);

A1 = [1];
A2 = [1];
b1 = [1];
b2 = [-1];
c1=[];
c2=[];
                
       
    flaghess = 'ifdois'; multhess = 1;
    Geral
    kvec = k;
    tvec = maint;
    evec = erro;
    
    kvecjacobi = kjacobi;
    tvecjacobi = tjacobi;
    evecjacobi = ejacobi;
    
    flaghess = 'jacobiintern'; x0=ones(n+m,1);
    Geral
    kvecintern = k;
    tvecintern = maint;
    evecintern = erro;
    
    flaghess = 'jacobidois'; x0=ones(n+m,1);
    Geral
    kvecdois = k;
    tvecdois = maint;
    evecdois = erro;
    
    x0=ones(n+m,1);
    [~,kvecyuan,tvecyuan,evecyuan] = Yuan_Varios_Fun(x0,n,m,A1,b1,c1,A2,b2,c2,mystr);
    
    asalvarstr = strcat('.\Salvar\resultsquadtestes1',num2str(n),num2str(m));
    save(asalvarstr,'kvec','tvec','evec','kvecjacobi','tvecjacobi','evecjacobi','kvecintern','tvecintern','evecintern','kvecdois','tvecdois','evecdois','kvecyuan','tvecyuan','evecyuan');



clear


alfa = 1/2; 
tol = 1e-4; kmax=1000;
minvaltol = 1;
maxvaltol=1e10;
%30 problemas de n=m=10
n=10; m=10;
mystr = 'quadtestes1';
x0 = ones(n+m,1);

A1 = [1];
A2 = [1];
b1 = [1];
b2 = [-1];
c1=[];
c2=[];
                
       
    flaghess = 'ifdois'; multhess = 1;
    Geral
    kvec = k;
    tvec = maint;
    evec = erro;
    
    kvecjacobi = kjacobi;
    tvecjacobi = tjacobi;
    evecjacobi = ejacobi;
    
    flaghess = 'jacobiintern'; x0=ones(n+m,1);
    Geral
    kvecintern = k;
    tvecintern = maint;
    evecintern = erro;
    
    flaghess = 'jacobidois'; x0=ones(n+m,1);
    Geral
    kvecdois = k;
    tvecdois = maint;
    evecdois = erro;
    
    x0=ones(n+m,1);
    [~,kvecyuan,tvecyuan,evecyuan] = Yuan_Varios_Fun(x0,n,m,A1,b1,c1,A2,b2,c2,mystr);
    
    asalvarstr = strcat('.\Salvar\resultsquadtestes1',num2str(n),num2str(m));
    save(asalvarstr,'kvec','tvec','evec','kvecjacobi','tvecjacobi','evecjacobi','kvecintern','tvecintern','evecintern','kvecdois','tvecdois','evecdois','kvecyuan','tvecyuan','evecyuan');



clear


alfa = 1/2; 
tol = 1e-4; kmax=1000;
minvaltol = 1;
maxvaltol=1e10;
%30 problemas de n=m=10
n=50; m=50;
mystr = 'quadtestes1';
x0 = ones(n+m,1);

A1 = [1];
A2 = [1];
b1 = [1];
b2 = [-1];
c1=[];
c2=[];
                
       
    flaghess = 'ifdois'; multhess = 1;
    Geral
    kvec = k;
    tvec = maint;
    evec = erro;
    
    kvecjacobi = kjacobi;
    tvecjacobi = tjacobi;
    evecjacobi = ejacobi;
    
    flaghess = 'jacobiintern'; x0=ones(n+m,1);
    Geral
    kvecintern = k;
    tvecintern = maint;
    evecintern = erro;
    
    flaghess = 'jacobidois'; x0=ones(n+m,1);
    Geral
    kvecdois = k;
    tvecdois = maint;
    evecdois = erro;
    
    x0=ones(n+m,1);
    [~,kvecyuan,tvecyuan,evecyuan] = Yuan_Varios_Fun(x0,n,m,A1,b1,c1,A2,b2,c2,mystr);
    
    asalvarstr = strcat('.\Salvar\resultsquadtestes1',num2str(n),num2str(m));
    save(asalvarstr,'kvec','tvec','evec','kvecjacobi','tvecjacobi','evecjacobi','kvecintern','tvecintern','evecintern','kvecdois','tvecdois','evecdois','kvecyuan','tvecyuan','evecyuan');



clear


alfa = 1/2; 
tol = 1e-4; kmax=1000;
minvaltol = 1;
maxvaltol=1e10;
%30 problemas de n=m=10
n=100; m=100;
mystr = 'quadtestes1';
x0 = ones(n+m,1);

A1 = [1];
A2 = [1];
b1 = [1];
b2 = [-1];
c1=[];
c2=[];
                
       
    flaghess = 'ifdois'; multhess = 1;
    Geral
    kvec = k;
    tvec = maint;
    evec = erro;
    
    kvecjacobi = kjacobi;
    tvecjacobi = tjacobi;
    evecjacobi = ejacobi;
    
    flaghess = 'jacobiintern'; x0=ones(n+m,1);
    Geral
    kvecintern = k;
    tvecintern = maint;
    evecintern = erro;
    
    flaghess = 'jacobidois'; x0=ones(n+m,1);
    Geral
    kvecdois = k;
    tvecdois = maint;
    evecdois = erro;
    
    x0=ones(n+m,1);
    [~,kvecyuan,tvecyuan,evecyuan] = Yuan_Varios_Fun(x0,n,m,A1,b1,c1,A2,b2,c2,mystr);
    
    asalvarstr = strcat('.\Salvar\resultsquadtestes1',num2str(n),num2str(m));
    save(asalvarstr,'kvec','tvec','evec','kvecjacobi','tvecjacobi','evecjacobi','kvecintern','tvecintern','evecintern','kvecdois','tvecdois','evecdois','kvecyuan','tvecyuan','evecyuan');



clear


alfa = 1/2; 
tol = 1e-4; kmax=1000;
minvaltol = 1;
maxvaltol=1e10;
%30 problemas de n=m=10
n=200; m=200;
mystr = 'quadtestes1';
x0 = ones(n+m,1);

A1 = [1];
A2 = [1];
b1 = [1];
b2 = [-1];
c1=[];
c2=[];
                
       
    flaghess = 'ifdois'; multhess = 1;
    Geral
    kvec = k;
    tvec = maint;
    evec = erro;
    
    kvecjacobi = kjacobi;
    tvecjacobi = tjacobi;
    evecjacobi = ejacobi;
    
    flaghess = 'jacobiintern'; x0=ones(n+m,1);
    Geral
    kvecintern = k;
    tvecintern = maint;
    evecintern = erro;
    
    flaghess = 'jacobidois'; x0=ones(n+m,1);
    Geral
    kvecdois = k;
    tvecdois = maint;
    evecdois = erro;
    
    x0=ones(n+m,1);
    [~,kvecyuan,tvecyuan,evecyuan] = Yuan_Varios_Fun(x0,n,m,A1,b1,c1,A2,b2,c2,mystr);
    
    asalvarstr = strcat('.\Salvar\resultsquadtestes1',num2str(n),num2str(m));
    save(asalvarstr,'kvec','tvec','evec','kvecjacobi','tvecjacobi','evecjacobi','kvecintern','tvecintern','evecintern','kvecdois','tvecdois','evecdois','kvecyuan','tvecyuan','evecyuan');
    


  
