-*

    This file contains functions needed to perform monodromy 
    (not performing monodromy itself).

*-

load 'parametrisation.m2';

-- Function takes dh parameters and returns {problem, solution} seed pair
findSeed = (alpha, r, d, theta) -> (
    -- Define local ring
    R := CC[x,y,z];

    -- Perform forward kinematics
    F := forwardKinematics(alpha, r, d, theta);
    ax := F_(0,3);
    ay := F_(1,3);
    az := F_(2,3);


    --Challenging to find seed in general, may want to take wrist partitioned approach! 

    S := submatrix'(F, {3}, {3});

    A := constructCayleySO();

    equations := {};
    for i from 0 to 2 do(
        for j from 0 to 2 do( 
            equations = append(equations, (x^2+y^2+z^2+1)*S_(i,j) - A_(i,j));
        );
    );

    Sols := solveSystem equations;

    R = CC[x, y, z, tx, ty, tz][ct1, ct2, ct3, ct4, ct5, ct6, st1, st2, st3, st4, st5, st6];

    problem := point {{(coordinates(Sols#0))#0, (coordinates(Sols#0))#1, ((coordinates(Sols#0))#2), ax, ay, az_CC}};
    solution := point {{cos(theta#0), cos(theta#1), cos(theta#2),cos(theta#3), cos(theta#4), cos(theta#5), sin(theta#0), sin(theta#1), sin(theta#2),sin(theta#3), sin(theta#4), sin(theta#5)_CC}};

    return {problem, solution};
);

-- Function returns list of 6 indexes of matrix A, s.t. A is 
--  non singular without those 6 rows. eps - presicion 
--  (if sigma_n > eps -> non singular -> good rows)
findBadRows = (A, eps) -> (
    -- 6 choose 12 subsets
    indexLst = subsets({0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11}, 6);
    
    -- Go through all possible indexes
    for i from 0 to #indexLst-1 do(
        indToDelete = indexLst#i;
        subMat = submatrix'(A, indToDelete, );  -- Delete rows
        decomposition = SVD(subMat);            -- Find SVD
        sigmas = decomposition#0;               -- Sigmas column
        if sigmas#-1 > eps then return indToDelete;  -- If non singular then return result
    );
    
   -- Return empty list if matrix is singular in any case
   return {};
);