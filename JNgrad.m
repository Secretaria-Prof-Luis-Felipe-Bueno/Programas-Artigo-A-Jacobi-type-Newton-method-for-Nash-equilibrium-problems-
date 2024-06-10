function g=JNgrad(x,i,teste)

if teste==1
    if i==1
        g=x{1}(1)-x{2}(1)^2;
    end
    
    if i==2
        g=x{2}(1)-x{1}(1)^2;
    end
end


if teste==2
    if i==1
        g=(x{1}(1)+10*x{2}(1)-21)*2;
    end
    
    if i==2
        g=(5*x{1}(1)+x{2}(1)-7)*2;
    end    
end    


if teste==3  
    if i==1
        g=-x{1}(1)^2*x{2}(1)^2+x{1}(1);
    end
    
    if i==2
        g=-x{2}(1)^2*x{1}(1)^2+x{2}(1);
    end    
end

if teste==4
    a=-x{1}(1)-x{2}(1)-x{3}(1);
    g=a+x{i}(1)+2*(max(x{i}(1)-1,0)-max(-x{i}(1)-1,0));
end
