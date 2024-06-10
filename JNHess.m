function H=JNHess(x,i,j,teste)

if teste==1
    if i==1
        if j==1
            H=1;
        else
            H=-2*x{2}(1);
        end
    end
    
    if i==2
        if j==1
            H=-2*x{1}(1);
        else
            H=1;
        end
    end
end


if teste==2
    if i==1
        if j==1
            H=2;
        else
            H=20;
        end
    end
    
    if i==2
        if j==1
            H=10;
        else
            H=2;
        end
    end
end




if teste==3
    if i==1
        if j==1
            H=-x{1}(1)*2*x{2}(1)^2+1;
        else
            H=-x{1}(1)^2*x{2}(1)*2;
        end
    end
    
    if i==2
        if j==1
            H=-x{2}(1)^2*x{1}(1)*2;
        else
            H=-x{2}(1)*2*x{1}(1)^2+1;
        end
    end
end

if teste==4
    if i==j
        if (x{i}(1)>-1 && x{i}(1)<1)
            H=0;
        else
            H=2;
        end
    else
        H=-1;
    end
end


  