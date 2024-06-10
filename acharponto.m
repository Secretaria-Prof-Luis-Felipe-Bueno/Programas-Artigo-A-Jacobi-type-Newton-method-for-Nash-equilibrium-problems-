myflag=0;

for i=-2:0.25:2
 for j=-2:0.25:2
  %for k=-2:0.25:2
   %for l=-2:0.25:2
    if (myflag==0)
     %c=avaliaden([i;j],[k;l]);
     c=avaliaden([i],[j]);
     %if (sum(sum(c>=0))==6)&&(sum(sum(c>0))>=1)
     if (sum(sum(c>=0))==2)&&(sum(sum(c>0))>=1)
        myflag=1;
        xgo=1;
        ygo=j;
        %xgo=[i;j];
        %ygo=[k,l];
     end
    end
   %end
  %end
 end
end

