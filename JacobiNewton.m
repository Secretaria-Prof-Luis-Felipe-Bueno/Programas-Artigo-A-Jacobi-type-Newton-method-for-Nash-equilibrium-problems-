function [x,it]=JacobiNewton(x,teste)

itmax=10;
tol=10^-6;
alpha=10^-1;
n=size(x,1);
ni=zeros(n,1);
g=cell(n,1);
H=cell(n,n);

it=0;
stop=0;
ntotal=0;
for i=1:n
    ni(i)=size(x{i},1);
    g{i}=JNgrad(x,i,teste);
    stop=stop+norm(g{i});
    ntotal=ntotal+ni(i);
end


M=zeros(ntotal);
Maux=zeros(ntotal);
b=zeros(ntotal,1);

while (stop>tol && it<itmax)
    
    M=[];
    b=[];
    for i=1:n
        Maux=[];
        H{i,i}=JNHess(x,i,i,teste);
        [~,p]=chol(H{i,i});
        if p>0
            H{i,i}=eye(ni(i));
        end
        
        for j=1:i-1
            H{i,j}=JNHess(x,i,j,teste);
            Maux=[Maux H{i,j}];
        end
        
        Maux=[Maux H{i,i}];
        
        for j=i+1:n
            H{i,j}=JNHess(x,i,j,teste);
            Maux=[Maux H{i,j}];
        end
        M=[M;Maux];
        
        b=[b;-g{i}];
    end
    
    
    M
    b
    [d,r]=linsolve(M,b)
    if r<0.001
        red=-1;
    else
        s=0;
        
        for i=1:n
            di{i}=d(s+1:s+ni(i));
            xT{i}=x{i}+di{i};
            s=s+ni(i);
        end
        
        red=1;
        
        for i=1:n
            fT{i}=JNfun(xT,i,teste);
            xaux=xT;
            xaux{i}=x{i};
            f{i}=JNfun(xaux,i,teste);
            red=min(red,f{i}-fT{i}-alpha*norm(di{i})^2);
        end
    end
    t=1;
    
    while red<0
        t=t/2;
        
        M=[];
        for i=1:n
            Maux=[];
            
            for j=1:i-1
                Maux=[Maux t*H{i,j}];
            end
            
            Maux=[Maux H{i,i}];
            
            for j=i+1:n
                Maux=[Maux t*H{i,j}];
            end
            M=[M;Maux];
        end
        
        M
        b
        [d,r]=linsolve(M,b)
        if r<0.001
            red=-1;
        else
            
            s=0;
            
            for i=1:n
                di{i}=d(s+1:s+ni(i));
                xT{i}=x{i}+t*di{i};
                s=s+ni(i);
            end
            
            red=1;
            
            for i=1:n
                fT{i}=JNfun(xT,i,teste);
                xaux=xT;
                xaux{i}=x{i};
                f{i}=JNfun(xaux,i,teste);
                red=min(red,f{i}-fT{i}-alpha*norm(di{i})^2);
            end
        end
    end
    
    x=xT;
    
    it=it+1;
    stop=0;
    for i=1:n
        g{i}=JNgrad(x,i,teste);
        stop=stop+norm(g{i});
    end
    x
    
end
