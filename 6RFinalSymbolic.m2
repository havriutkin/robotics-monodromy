
needsPackage "NumericalAlgebraicGeometry";
needsPackage "MonodromySolver";

 R = ZZ/7772777[x, y, z, tx, ty, tz, S1, S2, S3, S4, S5, S6, C1, C2, C3, C4, C5, C6, D1, D2, D3, D4, D5, D6, R1, R2, R3, R4, R5, R6][ct1, ct2, ct3, ct4, ct5, ct6, st1, st2, st3, st4, st5, st6]; -- we also need to add r's. 
-- Replaced ci si with the angles themselves. Parametrization is linear, but still we have problems.
-- Function to generate random DH parameters

constructEquations = () -> (
    ctheta := {ct1, ct2, ct3, ct4, ct5, ct6};
    stheta := {st1, st2, st3, st4, st5, st6};
    calpha := {C1, C2, C3, C4, C5, C6};
    salpha := {S1, S2, S3, S4, S5, S6};
    R := {R1, R2, R3, R4, R5, R6};
    D := {D1, D2, D3, D4, D5, D6};
   
    
    B := {tx, ty, tz};
    A := matrix{{-x^2-y^2+z^2+1, -2*y*z-2*x, 2*x*z-2*y}, {-2*y*z+2*x, -x^2+y^2-z^2+1, -2*x*y-2*z}, { 2*x*z+2*y, -2*x*y+2*z,x^2-y^2-z^2+1}};
    
    -- Construct SE(3) with unknowns
    T := i -> matrix{{ctheta#i, -stheta#i*calpha#i,  stheta#i*salpha#i, R#i*ctheta#i},
                    {stheta#i,  ctheta#i*calpha#i,  -ctheta#i*salpha#i, R#i*stheta#i},
                    {0,         salpha#i,            calpha#i,          D#i},
                    {0,         0,                      0,               1}};
    
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


end 
-*
restart

load "6RFinalSymbolic.m2"
needsPackage "MonodromySolver"
needsPackage "NumericalAlgebraicGeometry"
errorDepth=0;


DH = {{1,0,1,-1,1,0}, {0,1,0,0,0,1}, {0.675,0,0,0,0.6,0.140}, {0.3,0.65,0.155,0,0,0}}; --ENTER SPECIFIC DH PARAMETERS HERE IN THE FOLLOWING ORDER: SIN(ALPHA), COS(ALPHA), D, R
Params = {3.61498,2.35762e-15-5.95156e-15*ii,-1.9952,-0.652129,0.48607,0.445194}; --x, y, z, tx, ty, tz for specific problem
F = constructEquations();
-- Creating the ring map f
L = CC[ct1, ct2, ct3, ct4, ct5, ct6, st1, st2, st3, st4, st5, st6];
--{ct1=>ct1, ct2=>ct2, ct3=>ct3, ct4=>ct4, ct5=>ct5, ct6=>ct6, st1=>st1, st2=>st2, st3=>st3, st4=>st4, st5=>st5, st6=>st6, S1=>DH#0#0, S2=>DH#0#1, S3=>DH#0#2, S4=>DH#0#3,S5=>DH#0#4, S6=>DH#0#5, C1=>DH#1#0, C2=>DH#1#1, C3=>DH#1#2, C4=>DH#1#3,C5=>DH#1#4, C6=>DH#1#5, D1=>DH#2#0, D2=>DH#2#1, D3=>DH#2#2, D4=>DH#2#3,D5=>DH#2#4,D6=>DH#2#5,R1=>DH#3#0, R2=>DH#3#1, R3=>DH#3#2,R4=>DH#3#3,R5=>DH#3#4,R6=>DH#3#5, x=>Params#0, y=>Params#1, z=>Params#2, tx=>Params#3, ty=>Params#4, tz=>Params#5}
unknowns = {ct1, ct2, ct3, ct4, ct5, ct6, st1, st2, st3, st4, st5, st6};
specification = flatten{unknowns,Params,flatten(DH)}
f = map(L,R,specification);
-- Apply the ring map
func = (poly) -> (f(poly));
lst = apply(F,func)

I=ideal(lst);
dim(I)
degree(I)
*-
