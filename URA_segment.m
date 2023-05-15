function URAseg= URA_segment(D4th)
% To find the URA segment, you can make a matrix of ones that has a size
% equal to the aperture of the 4th order coarray, ,ultiply it by the fourth
% order coarray and in the holes you will get zeros. Of course, matlab
% realization will be different from this, but consider this concept. After
% you get the holes, calculate the lowest x and y values and take one step
% down to get the dimensions of the URA



% Generate 4th order coarray
% find maximum x and y values
% generate a ura large segment that have the same aperture size of the 4th order
% coarray
% find the elements present in the large ura segment but not in 4th order coarray
% find the minimum values of x and and y that is not common between the two arrays
% and generate a ura segment of aperture size (-x:x) * (-y:y)
%you need to find the distance of the closest point to the origin and based
%on that determine your ura segment

xapp=max(abs(D4th(:,1)));
yapp=max(abs(D4th(:,2)));
%1- Generate the URA segment
Max_URA=small_ura(xapp+1,yapp+1);
%determine the holes in the coarray
holes= setdiff(Max_URA,D4th,'rows');
distances= sqrt(holes(:,1).^2 +holes(:,2).^2);
URAseg= Max_URA;
[~,closest]= min(distances);
 xhole= floor(distances(closest));
 yhole= floor(distances(closest));   



while ~isempty(holes)
    distances= sqrt(holes(:,1).^2 +holes(:,2).^2);
    [~,closest]= min(distances);
    
    switch abs(holes(closest,1))>abs(holes(closest,2))
        case 1
            xhole= abs(holes(closest,1))-1;
            
        case 0
            yhole= abs(holes(closest,2))-1;
            
    end
    
    URAseg=small_ura(xhole,yhole);
    
    holes= intersect(URAseg,holes,'rows');
    
end


end