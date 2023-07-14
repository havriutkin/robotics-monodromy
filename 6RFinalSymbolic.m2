
needsPackage "NumericalAlgebraicGeometry";
needsPackage "MonodromySolver";

 R = QQ[x, y, z, tx, ty, tz, S1, S2, S3, S4, S5, S6, C1, C2, C3, C4, C5, C6, D1, D2, D3, D4, D5, D6, R1, R2, R3, R4, R5, R6][ct1, ct2, ct3, ct4, ct5, ct6, st1, st2, st3, st4, st5, st6]; -- we also need to add r's. 
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

restart

load "6RFinalSymbolic.m2"
needsPackage "MonodromySolver"
needsPackage "NumericalAlgebraicGeometry"
errorDepth=0;


DH = {{3,3,3,3,3,3}, {1,1,1,1,1,1}, {2,2,2,2,2,2}, {1,1,1,1,1,1}}; --ENTER SPECIFIC DH PARAMETERS HERE IN THE FOLLOWING ORDER: SIN(ALPHA), COS(ALPHA), D, R
Params = {1/8,1/8,1/8,1/8,1/8,1/8}; --x, y, z, tx, ty, tz for specific problem
F = constructEquations();
-- Creating the ring map f
L = QQ[ct1, ct2, ct3, ct4, ct5, ct6, st1, st2, st3, st4, st5, st6];
--{ct1=>ct1, ct2=>ct2, ct3=>ct3, ct4=>ct4, ct5=>ct5, ct6=>ct6, st1=>st1, st2=>st2, st3=>st3, st4=>st4, st5=>st5, st6=>st6, S1=>DH#0#0, S2=>DH#0#1, S3=>DH#0#2, S4=>DH#0#3,S5=>DH#0#4, S6=>DH#0#5, C1=>DH#1#0, C2=>DH#1#1, C3=>DH#1#2, C4=>DH#1#3,C5=>DH#1#4, C6=>DH#1#5, D1=>DH#2#0, D2=>DH#2#1, D3=>DH#2#2, D4=>DH#2#3,D5=>DH#2#4,D6=>DH#2#5,R1=>DH#3#0, R2=>DH#3#1, R3=>DH#3#2,R4=>DH#3#3,R5=>DH#3#4,R6=>DH#3#5, x=>Params#0, y=>Params#1, z=>Params#2, tx=>Params#3, ty=>Params#4, tz=>Params#5}
unknowns = {ct1, ct2, ct3, ct4, ct5, ct6, st1, st2, st3, st4, st5, st6};
specification = flatten{unknowns,Params,flatten(DH)}
f = map(L,R,specification);
-- Apply the ring map
func = (poly) -> (f(poly));
lst = apply(F,func)

-*
F1= f(F#0);
F2= f(F#1);
F3= f(F#2);
F4= f(F#3);
F5= f(F#4);
F6= f(F#5);
F7= f(F#6);
F8= f(F#7);
F9= f(F#8);
F10= f(F#9);
F11= f(F#10);
F12= f(F#11);
F13= f(F#12);
F14= f(F#13);
F15= f(F#14);
F16= f(F#15);
F17= f(F#16);
F18= f(F#17);

R=QQ[ct1,ct2,ct3,ct4,ct5,ct6,st1,st2,st3,st4,st5,st6];


S={F1,F2,F3,F4,F5,F6,F7,F8,F9,F10,F11,F12,F13,F14,F15,F16,F17};
*-
I=ideal(lst);
dim(I)
degree(I)
