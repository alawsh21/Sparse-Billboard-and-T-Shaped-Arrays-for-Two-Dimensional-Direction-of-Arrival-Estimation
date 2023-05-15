function d= Difference_Coarray(s)
lens= length(s(:,1));  
d(lens^2,2)= 0;


%Generate the second order coarray
for i= 1:lens
    index= i-1;    
    d((1:lens)+index*lens,:)= s(i,:) - s ; 
          
%     
%     for j= 1: lens
%        d(index*lens + j,:)= s(i,:) - s(j,:) ; 
%           
%     end 
        
end








end 




