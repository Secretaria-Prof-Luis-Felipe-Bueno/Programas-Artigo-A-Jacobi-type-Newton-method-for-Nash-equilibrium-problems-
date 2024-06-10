function [x,lambda,clak,erro1,erro2] = claudia(x0,n,m,A1,b1,c1,mystr,r1,r2)

%pegamos metade pq ? so a 1a parte n?
x0 = x0(1:2);

lambda = 0;
cx = x0(1) + sqrt(1+x0(2)^2);
gx = 2*A1*x0 + b1;
ax = [1, x0(2)/sqrt(1+x0(2)^2)];

clak=0; aux=0; kmax=1000;
gammamax = 1; gamma = gammamax;

while(aux==0)
    if clak>=kmax
        aux=1;
    end
        
    erro1 = norm(cx);
    gradlam = gx + lambda*ax';
    erro2 = norm(gradlam);
    if (erro1<1e-3)&&(erro2<1e-3)
        aux = 1;
    else
    
    Lk=2*A1; Lk(2,2) = Lk(2,2)+lambda/((1+x0(2)^2)^(3/2));
    bigmatrix = [Lk ax'; ax 0];
    biggrad = [-gradlam; -cx];
    vec = bigmatrix\biggrad;
    d = vec(1:2);
    lambda=vec(3);
    
    %update do sigma
    if gamma>=1.1*(lambda+gammamax)
        gamma = (gamma + gammamax + lambda)/2;
    else
        if gamma < lambda + gammamax
            gamma = max(1.5*gamma,lambda+gammamax);
        end
        
    end
    
    
    %update alfa
    condalfa=0; alfa=1;
    while condalfa==0
        diferenca = (x0+alfa*d)'*(A1*(x0+alfa*d)) + b1'*(alfa*d) - x0'*(A1*x0) + gamma*(norm([x0(1)+d(1) + sqrt(1+(x0(2)+d(2))^2)]) - norm(cx));
        if diferenca <=0
            condalfa=1;
        else
           alfa=alfa/2; 
        end
    end
    
    x0 = x0+alfa*d;
    
    cx = x0(1) + sqrt(1+x0(2)^2);
    gx = 2*A1*x0 + b1;
    ax = [1, x0(2)/sqrt(1+x0(2)^2)];
    
    clak = clak+1;
    
    end
    
end

x=x0;

end

