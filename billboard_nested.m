function S= billboard_nested(N1,N2,rotation)

d2 = [ 0 : N2-1] * (N1+1);
d1 = [ 0 : N1-1]+1+ d2(end) ;
d3 = [ 0 : N1-1] ;
d4 = [ 0 : N2-1] * (N1+1)+1+d1(end);

pos= unique([d1,d2]'); %rotated
pos2= unique([d3,d4]'); %not rotated


if rotation
S= [pos,zeros(length(pos2),1);...
    pos,pos;...
    zeros(length(pos2),1),pos];

else
   S= [pos2,zeros(length(pos2),1);...
       pos2,pos2;...
    zeros(length(pos2),1),pos2];
 
    
end 

S= unique(S,'rows');


end 