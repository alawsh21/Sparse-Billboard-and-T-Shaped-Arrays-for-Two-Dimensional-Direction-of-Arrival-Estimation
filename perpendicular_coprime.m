function S= perpendicular_coprime(M,N)
d1= [0:N-1]*M;
d2= [0:M-1]*N;

pos= unique([d1,d2]');

S= [pos,zeros(length(pos),1);...
    -pos,zeros(length(pos),1);...
    zeros(length(pos),1),pos];
    
S= unique(S,'rows');




end 