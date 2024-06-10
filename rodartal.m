t=[-0.2:0.01:0]; sv=t;
for i=1:length(t)
    sv(i) = Feval([t(i)-tx;x0(2)-tx],A1,b1,c1,0,'aplbuenovdd',1);
end
scatter(t,sv)