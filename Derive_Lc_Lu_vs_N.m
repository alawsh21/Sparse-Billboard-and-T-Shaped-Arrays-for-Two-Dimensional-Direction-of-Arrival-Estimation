clear 
close all
clc
%% Call Difference_Coarrays_Generation
Indicator=1;
Plot_Enable=1;

type_select=input(['1->B-Coprime,2->B-Nested,3->B-Super Nested,4->T-Coprime,5->T-Nested,6->T-Super Nested, type_select = ']);
select_operation=input(['0-> range of elements, 1-> user define, 2-> random test select_operation = ' ])


switch select_operation
    case 0 % Pre-defined range
        if type_select == 1 || type_select == 4 
%             N1=[2 3 3 4 3 5 5 6 5 7 7 8 7   9  9  10 9  11 11 12 11 13 13 14 13 15 15 16 16 15 17 17 18 17 19 19 19];
%             N2=[3 4 5 5 7 6 7 7 9 8 9 9 11  10 11 11 13 12 13 13 15 14 15 15 17 16 17 17 17 19 18 19 19 21 20 21 22];
            
            N1=[2 3 3 4 3 5 5 6 5 7 7 8 7   9  9];
            N2=[3 4 5 5 7 6 7 7 9 8 9 9 11  10 11];
            Nc = N1 + N2 - 1;
            N = 3*Nc - 2;
            disp('The total of coprime is, Nc = N1 + N2 -1')
            disp(['Nc = ' num2str(Nc) ' elements']) 
            disp('The total of 2D-coprime is, N = 3Nc - 2')   
            disp(['N = ' num2str(N) ' elements']) 
            N_Loop=length(Nc); 

        elseif type_select == 2 || type_select == 5
%             Nn = 4:40;
            Nn = 4:19;
            N_Loop=length(Nn);
            N = 3*Nn - 2;
            disp('The total of nested is, Nn = N1 + N2')
            disp(['Nn = ' num2str(Nn) ' elements']) 
            disp('The total of 2D-nested is, N = 3Nn - 2')   
            disp(['N = ' num2str(N) ' elements']) 
        elseif type_select == 3 || type_select == 6
%             Ns = 8:40;
            Ns = 8:19;
            Nn = Ns;
            N_Loop=length(Ns);
            N = 3*Ns - 2;
            disp('The total of super nested is, Ns = N1 + N2')
            disp(['Ns = ' num2str(Ns) ' elements']) 
            disp('The total of 2D-super nested is, N = 3Nn - 2')   
            disp(['N = ' num2str(N) ' elements']) 
        end
        
        
    case 1 % User Specified
        if type_select == 1 || type_select == 4  
            disp('Enter pairwise coprime integers, N1 < N2')
            N1 = input(['Enter N1, N1 = '])
            N2 = input(['Enter N2, N2 = '])
            disp('The total of coprime is, Nc = N1 + N2 -1')
            Nc = N1 + N2 - 1;
            N = 3*Nc - 2;
            disp(['Nc = ' num2str(Nc) ' elements']) 
            disp('The total of 2D-coprime is, N = 3Nc - 2')   
            disp(['N = ' num2str(N) ' elements']) 
            N_Loop=length(Nc);
        elseif type_select == 2 || type_select == 3 || type_select == 5 || type_select == 6
            disp('Enter the total of nested or super nested is, Nn = Ns = N1 + N2')
            disp('FOR SUPER NESTED, Ns >= 8 elements')
            disp('The total of 2D-Nested/Supernested is, N = 3Nn - 2 = 3Ns - 2')
            Nn = input(['Enter Nn, Nn = ']);
            N = 3*Nn - 2;
            disp('The total of 2D-Nested/Supernested is, N = 3Nn - 2')   
            disp(['N = ' num2str(N) ' elements']) 
            N_Loop=length(Nn);

        end
    case 2 %Random
        if type_select == 1 || type_select == 4  
            Nc = 2:22;
            Max_Num_Test = 10;
            Nc=round(Nc(1)+(Nc(end)-Nc(1))*rand(2,Max_Num_Test));   
            [~,ind_Coprime]= find(gcd(Nc(1,:),Nc(2,:)) == 1);            
            Nc=sort(Nc(:,ind_Coprime));
            [~,SUM_Nc]=sort(sum(Nc));
            Nc=Nc(:,SUM_Nc);
            N1=Nc(1,:);N2=Nc(2,:);
            N = 3*(N1+N2-1) - 2;
            disp('The total of 2D-Nested/Supernested is, N = 3Nc - 2') 
            disp(['N1 = ' num2str(N1) ' elements']) 
            disp(['N2 = ' num2str(N2) ' elements']) 
            disp(['N = ' num2str(N) ' elements']) 
            N_Loop=length(Nc);          
            
        elseif type_select == 2 || type_select == 5
%             Nn = 4:40;
            Nn = 4:19;
            Max_Num_Test = 10;
            Nn=unique(round(Nn(1)+(Nn(end)-Nn(1))*rand(1,Max_Num_Test)));              
            N = 3*Nn - 2;
            disp('The total of 2D-Nested/Supernested is, N = 3Nn - 2')   
            disp(['N = ' num2str(N) ' elements']) 
            N_Loop=length(Nn);
        elseif type_select == 3 || type_select == 6
%             Ns = 8:40;
            Ns = 8:19;
            Nn = Ns;            
            Max_Num_Test = 10;
            Nn=unique(round(Nn(1)+(Nn(end)-Nn(1))*rand(1,Max_Num_Test)));            
            N = 3*Nn - 2;
            disp('The total of 2D-Nested/Supernested is, N = 3Nn - 2')   
            disp(['N = ' num2str(N) ' elements']) 
            N_Loop=length(Nn);          
            
        end

        
end

if type_select == 2 || type_select == 3 || type_select == 5 || type_select == 6
    for ij=1:N_Loop
        if rem(Nn(ij),2) == 0
            N1(ij)= Nn(ij)/2;
            N2(ij)= Nn(ij)/2;
        else
            N1(ij)= (Nn(ij)-1)/2;
            N2(ij)= (Nn(ij)+1)/2;
        end
    end
end



TYPES = 6;
name = cell(TYPES, 1);

for type =  [type_select ]  % 1 : TYPES

    for ii=1:N_Loop
                   
        
        %% ========== Array geometries ==========
        switch (type)
            case 1 % Billboard
                name{type} = 'Billboard Coprime';
                S= billboard_coprime(N1(ii),N2(ii));                
                L = 3*(N1(ii)+N2(ii)-1)-2;
            case 2
                name{type} = 'Billboard Nested';
                S= billboard_nested(N1(ii),N2(ii),1);                
                L = 3*(N1(ii)+N2(ii))-2;        
            case 3
                name{type} = 'Billboard Super Nested';
                S= billboard_supernested(N1(ii),N2(ii));                
                L = 3*(N1(ii)+N2(ii))-2;  
            case 4 % T-Shaped
                name{type} = 'T-Shaped Coprime';            
                L = 3*(N1(ii)+N2(ii)-1)-2;
                S= perpendicular_coprime(N1(ii),N2(ii));             
            case 5
                name{type} = 'T-Shaped Nested';
                L = 3*(N1(ii)+N2(ii))-2;        
                S= perpendicular_nested(N1(ii),N2(ii),1);             
            case 6
                name{type} = 'T-Shaped Super Nested';
                L = 3*(N1(ii)+N2(ii))-2;        
                S= perpendicular_supernested(N1(ii),N2(ii)); 
            otherwise
        end
        disp(name{type})
        % ===================================== %

        %% Find the virtual lags through Simulation %%
        number_of_elements(ii)= length(S(:,1));
        D= unique(Difference_Coarray(S),'rows');    % Use S to find SODC
        D4th=unique(Difference_Coarray(D),'rows');  % Use SODC to find FODC
        URAseg= URA_segment(D4th);
        lc(ii)= length(URAseg(:,1)); % # Unique lags       
        lu(ii)=length(D4th(:,1)); % # Consecutive lags (uDOF)
        Num_Sensors(ii)=L;
        
        
        % [Second_Order_Diff_Co,Fourth_Order_Diff_Co,Holes]=Difference_Coarrays_Generation(S,Indicator, Plot_Enable);
        %% Find the virtual lags through Derivation %%
%%%%%%%%%%%%%%%%%%%%% Derivation of lc for all Arrays %%%%%%%%%%%%%%%%%%%%%
 % ====================================================================== %
           
        clear Edge

        switch type
 % ====================================================================== %
            case 1 % Billboard Coprime
                if (2*N1(ii) > N2(ii)) & (N1(ii) > 2) %& (gcd(N1(ii),N2(ii))==1)
                    Edge = N1(ii)*(N2(ii)+1)-1;  
                elseif (2*N1(ii) < N2(ii))  | N1(ii) == 2 & N2(ii) == 3 %& (gcd(N1(ii),N2(ii))==1)
                    Edge = N1(ii)*N2(ii)-N1(ii)+N2(ii)-1;
                end
                lc_Derived_B_Cop(ii)=(2*Edge+1)^2;    
                lc_Sim_B_Cop(ii)=lc(ii);  
                lu_Sim_B_Cop(ii)=lu(ii);  
                D_Derived_B_Cop(ii)=N1(ii)*(N2(ii)-1)*N1(ii)*(N2(ii)-1);
 % ====================================================================== %
            case 2 % Billboard Nested        
                Edge = (N1(ii)+1)*(N2(ii)+1)-2;
                lc_Derived_B_Nes(ii)=(2*Edge+1)^2;
                lc_Sim_B_Nes(ii)=lc(ii);    
                lu_Sim_B_Nes(ii)=lu(ii);    
                D_Derived_B_Nes(ii)=( N1(ii)+(N2(ii)-1)*(N1(ii)+1) )*(N1(ii)+(N2(ii)-1)*(N1(ii)+1));
 % ====================================================================== %
            case 3 % Billboard Super Nested    
                %Combined conditions
                if (rem(N1(ii),2) ~= 0) & (rem(N2(ii),2) ~= 0) | (rem(N1(ii),2) ~= 0) & (rem(N2(ii),2) == 0) & (rem(N2(ii)/2,2) == 0) | (rem(N1(ii),2) ~= 0) & (rem(N2(ii),2) == 0) & (rem(N2(ii)/2,2) ~= 0) 
                    Edge=N1(ii)*N2(ii)+N2(ii);                    
                elseif (rem(N1(ii)/2,2) ~= 0) & (rem(N2(ii)/2,2) ~= 0) & (rem(N1(ii),2) == 0) & (rem(N2(ii),2) == 0)
                    Edge = N1(ii)*(N2(ii)+1)+N1(ii)-1; 
                elseif (rem(N1(ii),2) == 0) & (rem(N2(ii),2) == 0) & (rem(N1(ii)/2,2) == 0) & (rem(N2(ii)/2,2) == 0)
                    Edge = N1(ii)*(N2(ii)+1)+N1(ii)/2-1; 
                elseif (rem(N1(ii),2) == 0) & (rem(N2(ii),2) ~= 0) & (rem(N1(ii)/2,2) == 0)
                    Edge = N1(ii)*(N2(ii)+1)+N1(ii)/2; 
                elseif (rem(N1(ii),2) == 0) & (rem(N2(ii),2) ~= 0) & (rem(N1(ii)/2,2) ~= 0)
                    Edge = N1(ii)*(N2(ii)+1)+N1(ii);        
                end
                lc_Derived_B_Sup(ii)=(2*Edge+1)^2; 
                lc_Sim_B_Sup(ii)=lc(ii);    
                lu_Sim_B_Sup(ii)=lu(ii);   
                D_Derived_B_Sup(ii)=(N1(ii)+(N2(ii)-1)*(N1(ii)+1))*(N1(ii)+(N2(ii)-1)*(N1(ii)+1));
 % ====================================================================== %
            case 4 % T-shaped Coprime
                if (N1(ii) ~=2) & (gcd(N1(ii),N2(ii)) == 1)
                    Edge=(N1(ii)+1)*N2(ii)+N1(ii)-1;
                elseif (N1(ii)== 2) & (gcd(N1(ii),N2(ii)) == 1)
                    Edge= 3*N2(ii)-1;
                end
                lc_Derived_T_Cop(ii)=(2*Edge+1)^2; 
                lc_Sim_T_Cop(ii)=lc(ii); 
                lu_Sim_T_Cop(ii)=lu(ii);    
                D_Derived_T_Cop(ii)=2*N1(ii)*(N2(ii)-1)*N1(ii)*(N2(ii)-1);
 % ====================================================================== %
            case 5 % T-shaped Nested
                Edge=2*(N1(ii)*N2(ii)+N2(ii)-1);
                lc_Derived_T_Nes(ii)=(2*Edge+1)^2; 
                lc_Sim_T_Nes(ii)=lc(ii);  
                lu_Sim_T_Nes(ii)=lu(ii);    
                D_Derived_T_Nes(ii)=2*(N1(ii)+(N2(ii)-1)*(N1(ii)+1))*(N1(ii)+(N2(ii)-1)*(N1(ii)+1));
 % ====================================================================== %
            case 6 % T-shaped Super Nested
                if (mod(N1(ii),2) ~= 0) && (mod(N2(ii),2) ~= 0) | (rem(N1(ii),2) ~= 0) & (rem(N2(ii),2) == 0) 
                    Edge = (N1(ii)+1)*(N2(ii)+1)+N1(ii)-1;
                elseif (rem(N1(ii)/2,2) ~= 0) & (rem(N2(ii)/2,2) ~= 0) & (rem(N1(ii),2) == 0) & (rem(N2(ii),2) == 0) |  (rem(N1(ii),2) == 0) & (rem(N2(ii),2) ~= 0) & (rem(N1(ii)/2,2) == 0)
                    Edge=N1(ii)*N2(ii)+5/2*N1(ii);
                elseif (rem(N1(ii),2) == 0) & (rem(N2(ii),2) == 0) & (rem(N1(ii)/2,2) == 0) & (rem(N2(ii)/2,2) == 0) 
                    Edge=N1(ii)*N2(ii)+5/2*N1(ii)-1; 
                elseif (rem(N1(ii),2) == 0) & (rem(N2(ii),2) ~= 0) & (rem(N1(ii)/2,2) ~= 0)
                    Edge = N1(ii)*N2(ii)+2*N1(ii)+(N2(ii)-1)/2+1; 
                end
                lc_Derived_T_Sup(ii)=(2*Edge+1)^2; 
                lc_Sim_T_Sup(ii)=lc(ii); 
                lu_Sim_T_Sup(ii)=lu(ii); 
                D_Derived_T_Sup(ii)=2*(N1(ii)+(N2(ii)-1)*(N1(ii)+1))*(N1(ii)+(N2(ii)-1)*(N1(ii)+1));
        end
        lc_Derived(ii)=(2*Edge+1)^2;
 % ====================================================================== %

    end % End of elements loop


%%


figure(100)

plot(Num_Sensors,lc,'bo',Num_Sensors,lc_Derived,'r+')
hold on
xlabel('Number of physical sensors, \itN');
ylabel('\itl_c \rm or (uDOF)')
legend('Simulated','Derived')
title([num2str(name{type})])

end

%%


[lc lc_Derived]
[Edge (lc-1)/2]


Edge_Sim=(sqrt(lc)-1)/2;
[lc lu  lc_Derived]

% figure
scatterplot(D4th) % plot the last FODC
hold on

if length(lc) == 1

rectangle('Position',[-Edge_Sim -Edge_Sim (2*Edge_Sim) (2*Edge_Sim)],'EdgeColor','b')

rectangle('Position',[-Edge -Edge (2*Edge) (2*Edge)],'EdgeColor','r')

end