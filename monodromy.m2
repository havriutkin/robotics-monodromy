-*

    This file contains functions needed to perform monodromy 
    (not performing monodromy itself).

*-

load "parametrisation.m2";

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

-- Function takes list of equations and list of points that represent solutions
--  returns sublist of solutions, that contains actual solutions to given equations
filterSolutions = (solutions, equations) -> (
    R := CC[ct1, ct2, ct3, ct4, ct5, ct6, st1, st2, st3, st4, st5, st6];    -- Declare local ring
    result = [];

    for i from 0 to length(solutions) do(
        flag := 1;                              -- If flag == 1 then solution is valid
        solution = coordinates(solutions#i);    -- Transofrm point to array
        solDict = {ct1 => solution#0,           -- Create dictionary
                    ct2 => solution#1, 
                    ct3 => solution#2, 
                    ct4 => solution#3, 
                    ct5 => solution#4, 
                    ct6 => solution#5, 
                    st1 => solution#6, 
                    st2 => solution#7, 
                    st3 => solution#8, 
                    st4 => solution#9, 
                    st5 => solution#10, 
                    st6 => solution#11,}

        -- Substitute solution to each equation
        for j from 0 to length(equations) do(
            if sub(equations#j, solDict) != 0 then flag := 0;
        );

        -- Add solution to result
        if flag == 1 then result = append(result, solution);
    );

    return result;
);