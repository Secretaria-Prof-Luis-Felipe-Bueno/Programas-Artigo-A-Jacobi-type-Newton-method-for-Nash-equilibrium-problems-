function [A,AL] = cholmod(A,tau)
%A quadrada, acha LLt ou senao de A+E conforme o artigo
%A deve sre simetrica
%A vai ser a nova matrix LLt e AL eh a matriz simetrica com L na diag de
%baixo e A na de cima, estilo cholesky matlab


faseum=1;
y = max(max(abs(A)));
n=length(A(1,:));
j=1;
L=A;
pivoteamento=1;

if A==zeros(n,n)
    A = tau*eye(n,n);
else
    
    %na fase 1, A eh pd
    while ((j<=n)&&(faseum==1))
        %acha indice do max diagonal de A
        i=j;
        for ind = j:n
            if A(ind,ind)>A(i,i)
                i=ind;
            end
        end
        if (i~=j)&&(pivoteamento==1) %swap rows and columns i and j of A
            vazio=A(i,1:n);vazio2=A(1:n,i);
            A(i,1:n) = A(j,1:n);
            A(j,1:n) = vazio;
            A(1:n,i) = A(1:n,j);
            A(1:n,j) = vazio2;
        end
        %computando min j+1:n aii-aij/ajj
        compara=zeros(1,n);
        for ind=j+1:n
            compara(ind)=A(ind,ind)-(A(ind,j)^2/A(j,j));
        end
        if min(compara(j+1:n))<tau*y
            faseum=0; %ai matriz nao eh pd
        else
            %nesse caso faz a j-esima parte da fatorizacao
            A(j,j)=sqrt(A(j,j)); %L
            for i=j+1:n
                A(i,j)=A(i,j)/A(j,j); %L
                for k=j+1:i
                    A(i,k)=A(i,k)-A(i,j)*A(k,j); %L
                end
            end
            j=j+1;
        end
        
    end %end while da fase 1
    
    
    %fase 2, matriz nao spd
    if faseum==0
        k=j-1; %numero de its da fase anterior
        
        %lower gerschgorin bounds
        g=zeros(1,n);
        maxg = k+1;
        for i=k+1:n
            g(i) = A(i,i) - sum(abs(A(i,k+1:i-1))) - sum(abs(A(i+1:n,i)));
            if g(i)>g(maxg) %indice do max g(i):
                maxg=i;
            end
        end
        
        delta=0;
        for j=k+1:n-2
            %pivoteia no max lower gersch bound
            if (maxg~=j)&&(pivoteamento==1)
                vazio=A(maxg,1:n);vazio2=A(1:n,maxg);
                A(maxg,1:n) = A(j,1:n);
                A(j,1:n) = vazio;
                A(1:n,maxg) = A(1:n,j);
                A(1:n,j) = vazio2;
            end
            
            %calcula ejj e soma na diag
            normj = sum(abs(A(j+1,n,j)));
            delta = max([0,delta,-A(j,j)+max([normj,tau*y])]);
            if delta>0
                A(j,j)=A(j,j)+delta;
            end
            
            %update gerschoring bounds:
            if A(j,j)~=normj
                temp=1-(normj/A(j,j));
                for i=j+1:n
                    g(i) = g(i) + abs(A(i,j))*temp;
                end
            end
            
            %j-esima fatorizacao
            A(j,j)=sqrt(A(j,j));
            for i=j+1:n
                A(i,j)=A(i,j)/A(j,j);
                for k=j+1:i
                    A(i,k)=A(i,k)-A(i,j)*A(k,j);
                end
            end
            
        end %end for do j
        
        %ultima submatriz para n-1 e n
        [lambdas] = eigs(A(n-1:n,n-1:n));
        lambinha=min(lambdas); lambao = max(lambdas);
        del = max([0,delta,-lambinha +tau*max([y,1/((1-tau)*(lambao-lambinha))]) ]);
        if del>0
            A(n-1,n-1) = A(n-1,n-1)+del;
            A(n,n) = A(n,n)+del;
            delta=del;
        end
        A(n-1,n-1) = sqrt(A(n-1,n-1));
        A(n,n-1) = A(n,n-1)/A(n-1,n-1);
        A(n,n) = sqrt(A(n,n) - (A(n,n-1))^2);
        
    end %end fase 2
    
    AL=A;
    
    %refazendo a A
    B = zeros(n,n);
    for i=1:n
        for j=1:i
            B(i,j)=A(i,j);
        end
    end
    B = B*B';
    A=B;
    
end


end

