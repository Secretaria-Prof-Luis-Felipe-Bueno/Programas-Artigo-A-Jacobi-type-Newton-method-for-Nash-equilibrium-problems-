%testepaper
n = 1; m = 1;
for i=1:30
    strmat = num2str(i);
    mystr = strcat('.\Salvar\matmistas',num2str(n),num2str(m),strmat);
    A1 = rand(n,m);
    %A1=(A1')*A1;
    A2 = rand(n,m);
    %A2=(A2')*A2;
    b1 = rand(n,1);
    b2 = rand(m,1);
    c1=rand(1,1);
    c2=rand(1,1);
    save(mystr,'A1','A2','b1','b2','c1','c2');    
end
n = 10; m = 10;
for i=1:30
    strmat = num2str(i);
    mystr = strcat('.\Salvar\matmistas',num2str(n),num2str(m),strmat);
    A1 = rand(n,m);
    %A1=(A1')*A1;
    A2 = rand(n,m);
    %A2=(A2')*A2;
    b1 = rand(n,1);
    b2 = rand(m,1);
    c1=rand(1,1);
    c2=rand(1,1);
    save(mystr,'A1','A2','b1','b2','c1','c2');     
end
n = 100; m = 100;
for i=1:30
    strmat = num2str(i);
    mystr = strcat('.\Salvar\matmistas',num2str(n),num2str(m),strmat);
    A1 = rand(n,m);
    %A1=(A1')*A1;
    A2 = rand(n,m);
    %A2=(A2')*A2;
    b1 = rand(n,1);
    b2 = rand(m,1);
    c1=rand(1,1);
    c2=rand(1,1);
    save(mystr,'A1','A2','b1','b2','c1','c2');     
end
n = 1000; m = 100;
for i=1:30
    strmat = num2str(i);
    mystr = strcat('.\Salvar\matmistas',num2str(n),num2str(m),strmat);
    A1 = rand(n,m);
    %A1=(A1')*A1;
    A2 = rand(n,m);
    %A2=(A2')*A2;
    b1 = rand(n,1);
    b2 = rand(m,1);
    c1=rand(1,1);
    c2=rand(1,1);
    save(mystr,'A1','A2','b1','b2','c1','c2');     
end
n = 1000; m = 1000;
for i=1:30
    strmat = num2str(i);
    mystr = strcat('.\Salvar\matmistas',num2str(n),num2str(m),strmat);
    A1 = rand(n,m);
    %A1=(A1')*A1;
    A2 = rand(n,m);
    %A2=(A2')*A2;
    b1 = rand(n,1);
    b2 = rand(m,1);
    c1=rand(1,1);
    c2=rand(1,1);
    save(mystr,'A1','A2','b1','b2','c1','c2');     
end
n = 500; m = 500;
for i=1:30
    strmat = num2str(i);
    mystr = strcat('.\Salvar\matmistas',num2str(n),num2str(m),strmat);
    A1 = rand(n,m);
    %A1=(A1')*A1;
    A2 = rand(n,m);
    %A2=(A2')*A2;
    b1 = rand(n,1);
    b2 = rand(m,1);
    c1=rand(1,1);
    c2=rand(1,1);
    save(mystr,'A1','A2','b1','b2','c1','c2');     
end
n = 50; m = 50;
for i=1:30
    strmat = num2str(i);
    mystr = strcat('.\Salvar\matmistas',num2str(n),num2str(m),strmat);
    A1 = rand(n,m);
    %A1=(A1')*A1;
    A2 = rand(n,m);
    %A2=(A2')*A2;
    b1 = rand(n,1);
    b2 = rand(m,1);
    c1=rand(1,1);
    c2=rand(1,1);
    save(mystr,'A1','A2','b1','b2','c1','c2');     
end