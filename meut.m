if strcmp(flaghess,'tfixo1')
        for i=1:length(xvec)
        scatvecy(i) = Feval([ty-x0(1);ty-yvec(i)],A1,b1,c1,0,'aplbuenovdd',1);
        end
        for i=1:100
        z(i) =  Feval([x0(1)-ty;t(i)-tx],A2,b2,c2,0,'aplbuenovdd',2);
        end
    end