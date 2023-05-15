function S= billboard_supernested(N1,N2)

pos = super_nested( N1, N2, 2 )-1;


S= [pos,zeros(length(pos),1);...
    zeros(length(pos),1),pos;pos,pos]; 
S= unique(S,'rows');


end 