function f=JNfun(x,i,teste)

if teste==1 %ptos eq. (0,0) e (1,1)
    if i==1
        f=0.5*x{1}(1)^2-x{2}(1)^2*x{1}(1);
    end
    
    if i==2
        f=0.5*x{2}(1)^2-x{1}(1)^2*x{2}(1);
    end    
end


if teste==2 %pto eq (1,2)
    if i==1
        f=(x{1}(1)+10*x{2}(1)-21)^2;
    end
    
    if i==2
        f=(5*x{1}(1)+x{2}(1)-7)^2;
    end    
end


if teste==3 %pto eq (0,0). Pto KKT não de eq (1,1) 
    if i==1
        f=-x{1}(1)^3*x{2}(1)^2/3+x{1}(1)^2/2;
    end
    
    if i==2
        f=-x{2}(1)^3*x{1}(1)^2/3+x{2}(1)^2/2;
    end    
end


if teste==4 %pto eq (0,0). Seis ptos de acumulação para Gauss Seidel, (1,1,-1), (1,-1,-1), (1,-1,1)... nenhum deles é pto eq. eq (1,1) 
    a= -x{1}(1)*x{2}(1)-x{1}(1)*x{3}(1)-x{2}(1)*x{2}(1);
    b= max(x{1}(1)-1,0)^2+max(x{2}(1)-1,0)^2+max(x{3}(1)-1,0)^2;
    c= max(-x{1}(1)-1,0)^2+max(-x{2}(1)-1,0)^2+max(-x{3}(1)-1,0)^2;
    f= a+b+c;
end