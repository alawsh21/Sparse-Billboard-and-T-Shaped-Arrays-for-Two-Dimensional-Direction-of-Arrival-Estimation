function S= billboard_coprime(M,N)


d1 = [ 0 : N-1] * M ;
d2 = [ 1 : M-1] * N;
pos = unique( sort([d1 d2]) ); % Coding by Dr. Jianyan LIU 
pos= pos';


S= [pos,zeros(length(pos),1);...
    zeros(length(pos),1),pos;pos,pos]; 
S= unique(S,'rows');

end
