function [NKKTJN, NKKTNP, NENJN, NENNP]=testealeatorio(n)

  NKKTJN=0;
  NKKTNP=0;
  NENJN=0;
  NENNP=0;

  tol=10^-4;
itmax=200;
hold on

for i=1:n
    x=-1+3*rand;
    y=-1+3*rand;
    
    [xJ,yJ,f1,f2,g1k,g2k,H1k,H2k]=testeJacobiNewton(x,y,2,0,itmax,tol);
    %plot(f1,f2,'*b')
    if g1k^2+g2k^2< tol
        NKKTJN=NKKTJN+1;
        if min(H1k,H2k) > 0
            NENJN=NENJN+1;
        end%if
    else
        x
        y
        g1k
        g2k
        pause
    end %if
    
    [xN,yN,f1,f2,g1k,g2k,H1k,H2k]=NewtonPuro(x,y,2,itmax,tol);
    %plot(f1,f2,'*r')
    if g1k^2+g2k^2< tol
        NKKTNP=NKKTNP+1;
        if min(H1k,H2k) > -tol
            if min(H1k,H2k) < tol
                xN
                yN
            end%if
            NENNP=NENNP+1;
        end%if
    end%if
end%for

end
