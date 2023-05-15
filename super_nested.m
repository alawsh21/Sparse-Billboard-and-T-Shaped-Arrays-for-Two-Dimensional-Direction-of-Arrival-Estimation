% Super-nested arrays given N1, N2, and Q
%
% Written by Chun-Lin Liu
% E-mail : cl.liu@caltech.edu
% Project website: http://systems.caltech.edu/dsp/students/clliu/SuperNested.html
%
% Our Papers:
% Journal
% [1] C.-L. Liu and P. P. Vaidyanathan, “Super Nested Arrays: Linear Sparse Arrays with Reduced Mutual Coupling - Part I: Fundamentals,” IEEE Transactions on Signal Processing, vol. 64, no. 15, pp. 3997-4012, Aug. 2016. 
% [2] C.-L. Liu and P. P. Vaidyanathan, “Super Nested Arrays: Linear Sparse Arrays with Reduced Mutual Coupling - Part II: High-Order Extensions,” IEEE Transactions on Signal Processing, vol. 64, no. 16, pp. 4203-4217, Aug. 2016. 
% Conference
% [3] C.-L. Liu and P. P. Vaidyanathan, “Super Nested Arrays: Sparse Arrays with Less Mutual Coupling than Nested Arrays,” in Proc. of 2016 IEEE International Conference on Acoustics Speech and Signal Processing (ICASSP 2016), pp. 2976-2980, Shanghai, China, Mar. 2016.
% [4] C.-L. Liu and P. P. Vaidyanathan, “High Order Super Nested Arrays,” in Proc. of the Ninth IEEE Sensor Array and Multichannel Signal Processing Workshop (SAM 2016), Rio de Janeiro, Brazil, Jul. 2016.
%
% Last revised on August 14, 2016
%


function S = super_nested( N1, N2, Q )
%SUPER_NESTED Returns sensor locations of super nested arrays
%   Given parameters N1, N2, and Q, return a super-nested array S (a column
%   vector of length N1+N2)

    if (Q == 1)
        % >>>>>>>>>>>>>>>>>>>> Nested arrays <<<<<<<<<<<<<<<<<<<<
        % Please see [5].
        % [5] P. Pal and P. P. Vaidyanathan, “Nested arrays: A novel approach to array processing with enhanced degrees of freedom,” IEEE Trans. Signal Process., vol. 58, no. 8, pp. 4167–4181, Aug 2010.
        S = [	1 : N1, ... % Dense ULA
				(1:N2)*(1+N1)].'; % Sparse ULA
    elseif (Q == 2)
        % >>>>>>>>>>>>>>>>>>>> Second-order super-nested arrays <<<<<<<<<<<<<<<<<<<<
        % Please see [1,3].
        % Generate (A1, B1, A2, B2) according to N1
        switch (mod(N1, 4))
            case 0
                r = N1 / 4;
                A1 = r; B1 = r - 1; A2 = r - 1; B2 = r - 2;
            case 1
                r = (N1 - 1) / 4;
                A1 = r; B1 = r - 1; A2 = r - 1; B2 = r - 1;
            case 2
                r = (N1 - 2) / 4;
                A1 = r + 1; B1 = r - 1; A2 = r; B2 = r - 2;
            case 3
                r = (N1 - 3)/4;
                A1 = r; B1 = r; A2 = r; B2 = r - 1;
            otherwise
                error('Error');
        end
        S = [   1 + 2 * (0 : A1), ...                   % X_1^{(Q)}
                (N1 + 1) - (1 + 2 * (0 : B1)), ...      % Y_1^{(Q)}
                (N1 + 1) + (2 + 2 * (0 : A2)), ...      % X_2^{(Q)}
                2 * (N1 + 1) - (2 + 2 * (0 : B2)), ...  % Y_2^{(Q)}
                (N1 + 1) * (2 : N2), ...                % Z_1^{(Q)}
                N2 * (N1 + 1) - 1                       % Z_2^{(Q)}
             ];
        S = sort(unique(S)).';
    else
        % >>>>>>>>>>>>>>>>>>>> Qth-order super-nested arrays <<<<<<<<<<<<<<<<<<<<
        % Please see [2,4].
        if ( mod(N1, 2) == 0)
            % N1 is EVEN
            S = super_nested( N1, N2, 2 );
            % Determine the cutoff between X_1^{(Q)} and Y_1^{(Q)}
            for nn = 1 : length(S)-1
                if (S(nn+1) - S(nn) == 1)
                    X = S(nn);
                    break;
                end
            end
            for q = 3 : Q
                % Select X_q^{(q)} and Y_q^{(q)}
                index_XQ = (S > (q-2)*(N1+1) & S <= (q-2)*(N1+1) + X);
                index_YQ = (S > (q-2)*(N1+1) + X & S < (q-1)*(N1+1));
                XQ = S( index_XQ );
                YQ = S( index_YQ );
                if (length(XQ) >= 3 && length(YQ) >= 3)
                    if ( mod(length(XQ), 2) == 0 )
                        % With extra term
                        XQ( 2 : 2 : end-2 ) = XQ( 2 : 2 : end-2 ) + (N1+1);
                    else
                        % No extra term
                        XQ( 2 : 2 : end-1 ) = XQ( 2 : 2 : end-1 ) + (N1+1);
                    end
                    
                    if ( mod(length(YQ), 2) == 0 )
                        % With extra term
                        YQ( end-1 : -2 : 3 ) = YQ( end-1 : -2 : 3 ) + (N1+1);
                    else
                        % No extra term
                        YQ( end-1 : -2 : 2 ) = YQ( end-1 : -2 : 2 ) + (N1+1);
                    end
                    S( S == (q-1)*(N1+1) ) = (N2 + 1 - (q - 1)) * (N1 + 1) - 2^(q-1) + 1;
                    S( index_XQ ) = XQ;
                    S( index_YQ ) = YQ;
                else
                    disp(['Only up to order-', num2str(q-1), ' is defined. Return the ', num2str(q-1), 'th-order super nested array'])
                    break;
                end
            end
            
            % Sort S
            S = sort(unique(S));
        else
            % N1 is ODD
            % Check Q satisfies our conditions or not
            if (N1 < 3 * 2^Q - 1 || N2 < 3 * Q - 4)
                disp('Warning! The resultant super-nested array might not be a restricted array!')
            end
            S = [];
            for q = 1 : Q - 1
                ell = 0 : floor(((N1 + 1) / 2^q - 1) / 2);
                S = [S,	(q - 1) * (N1 + 1) + 2^(q - 1) + 2^q * ell, q * (N1 + 1) - 2^(q - 1) - 2^q * ell, ...
                        (N2 + 1 - q) * (N1 + 1) - 2^q + 1];
            end
            ell = 0 : floor( (N1 + 1) / 2^Q - 1 );
            S = [S, (Q - 1) * (N1 + 1) + 2^(Q - 1) + 2^(Q - 1) * ell, Q * (N1 + 1) - 2^(Q - 1) - 2^(Q - 1) * ell];
            S = [S, (N1 + 1) * (Q : N2)];
            S = sort(unique(S)).';
        end
    end

end

