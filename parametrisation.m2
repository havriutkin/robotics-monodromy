-*

    This package contains functions related to parametrisation.

    Functions in this package: 
        - constructParametriseSE ():
            returns parametrised SE(3) matrix, parametrising rotation along each axies 
            using cx, sx, cy, sy, cz, sz. Translation parametrised by ax, ay, az.

        - constructRationalParametriseSE ():
            returns parametrised SE(3) matrix, parametrising rotation along each axies 
            using rational parametrisation (tx, ty, tz). Translation parametrised by ax, ay, az.

        - constructCayleySO ():
            returns parametrised SO(3) using Cayley parametrisation 
            multiplied by (x^2 + y^2 + z^2 + 1).

        - constructParametrisedEquations ():
            returns list of symbolical inverse kinematics equations 
            (dh parameters are symbolical)

*-

-- Function returns parametrised SE(3) using c_i and s_i
constructParametriseSE = () -> (
    -- Define local ring
    R := CC[ax, ay, az, cx, cy, cz, sx, sy, sz];

    -- SE(3) along X
    transX := matrix{{1, 0,  0,   ax},
                    {0, cx, -sx, 0},
                    {0, sx, cx,  0},
                    {0, 0,  0,   1}};
    
    -- SE(3) along Y
    transY := matrix{{cy,  0, sy, 0},
                    {0,   1, 0,  ay},
                    {-sy, 0, cy, 0},
                    {0,   0, 0,  1}};
                
    -- SE(3) along Z            
    transZ := matrix{{cz, -sz, 0, 0},
                    {sz, cz,  0, 0},
                    {0,  0,   1, az},
                    {0,  0,   0, 1}};
    
    result := transX * transY * transZ;
    return result;
);

-- Function returns parametrised SE(3) using rational parametrisation of cos and sin
constructRationalParametriseSE = () -> (
    -- Define local ring
    R := CC[ax, ay, az, tx, ty, tz];

    -- SE(3) along X
    transX := matrix{{1, 0, 0, ax}, 
                    {0, (1-tx^2)/(1+tx^2), (-2*tx)/(1+tx^2)}, 
                    {0, (2*tx)/(1+tx^2), (1-tx^2)/(1+tx^2), 0}, 
                    {0, 0, 0, 1}};
    
    -- SE(3) along Y
    transY := matrix{{(1-ty^2)/(1+ty^2), 0, (2*ty)/(1+ty^2), 0}, 
                    {0, 1, 0, ay},
                    {(-2*ty)/(1+ty^2), 0, (1-ty^2)/(1+ty^2), 0}, 
                    {0, 0, 0, 1}};
    
    -- SE(3) along Z
    transZ := matrix{{(1-tz^2)/(1+tz^2), (-2*tz)/(1+tz^2), 0, 0}, 
                    {(2*tz)/(1+tz^2), (1-tz^2)/(1+tz^2), 0, 0},
                    {0, 0, 1, az}, 
                    {0, 0, 0, 1}};
                
    result := transX * transY * transZ;
    return result;
);

-- Function returns parametrised SO(3) using Cayley parametrisation 
-- multiplied by (x^2 + y^2 + z^2 + 1)
constructCayleySO = () -> (
    -- Define local ring
    R := CC[x, y, z];

    A := matrix{{-x^2-y^2+z^2+1, -2*y*z-2*x,      2*x*z-2*y}, 
                {-2*y*z+2*x,     -x^2+y^2-z^2+1, -2*x*y-2*z}, 
                { 2*x*z+2*y,     -2*x*y+2*z,      x^2-y^2-z^2+1}};

    return A;
)

-- Function returns list of symbolical inverse kinematics equations 
-- (dh parameters are symbolical)
constructParametrisedEquations = () -> (
    -- Define local ring
    L := QQ[tx, ty, tz, S1, S2, S3, S4, S5, S6, C1, C2, C3, C4, C5, C6, D1, D2, D3, D4, D5, D6, R1, R2, R3, R4, R5, R6][ct1, ct2, ct3, ct4, ct5, ct6, st1, st2, st3, st4, st5, st6];

    -- Define symbolical DH parameters
    ctheta := {ct1, ct2, ct3, ct4, ct5, ct6};
    stheta := {st1, st2, st3, st4, st5, st6};
    calpha := {C1, C2, C3, C4, C5, C6};
    salpha := {S1, S2, S3, S4, S5, S6};
    R := {R1, R2, R3, R4, R5, R6};
    D := {D1, D2, D3, D4, D5, D6};
   
    -- Define symbolical SE(3) by cayley parametrisation 
    B := {tx, ty, tz};
    A := constructCayleySO();   -- This matrix is already multiplied by (x^2 + y^2 + z^2 + 1)
    
    -- Construct SE(3) with unknowns
    T := i -> matrix{{ctheta#i, -stheta#i*calpha#i,   stheta#i*salpha#i, R#i*ctheta#i},
                     {stheta#i,  ctheta#i*calpha#i,  -ctheta#i*salpha#i, R#i*stheta#i},
                     {0,         salpha#i,            calpha#i,          D#i},
                     {0,         0,                   0,                 1}};
    
    Tf := product(#D, i -> T(i));
    
    
    -- Set correspond elements of A and T to be equal
    myEquations := {};
    for i from 0 to 2 do(
        for j from 0 to 2 do( 
            myEquations = append(myEquations, (x^2+y^2+z^2+1)*Tf_(i,j) - A_(i,j));
        );
    );

    for k from 0 to 2 do(
            myEquations = append(myEquations, Tf_(k,3) - B#k); 
    );

    for l from 0 to 5 do(
    myEquations = append(myEquations, (ctheta#l)^2 + (stheta#l)^2 -1);
    );
    
    return myEquations;
);