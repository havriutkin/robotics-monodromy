R = ZZ/7772777[x, y, z, tx, ty, tz, S1, S2, S3, C1, C2, C3, D1, D2, D3, R1, R2, R3][ct1, ct2, ct3, st1, st2, st3]; -- we also need to add r's.

constructEquations = () -> (
    ctheta := {ct1, ct2, ct3};
    stheta := {st1, st2, st3};
    calpha := {C1, C2, C3};
    salpha := {S1, S2, S3};
    R := {R1, R2, R3};
    D := {D1, D2, D3};
   
    
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

    for l from 0 to 2 do(
    myEquations = append(myEquations, (ctheta#l)^2 + (stheta#l)^2 -1);
    );
    
    return myEquations;

);


end

restart

load "KukaExampleThreeArm.m2"
load "forwardKinematics.m2"
needsPackage "MonodromySolver"
needsPackage "NumericalAlgebraicGeometry"

S = CC[x,y,z];
--DH parameters for Kuka
alpha = {pi/2,0,pi/2};
d = {675,0,0};
r = {0,0,0};
theta = {pi/2,-pi/2,0};
-- SE(3) matrix for Kuka
M = forwardKinematics(alpha,r,d,theta)
-- x,y,z,tx,ty,tz for Kuka
txAns = sub(0,ZZ/7772777)
tyAns = sub(0,ZZ/7772777)
tzAns = sub(675,ZZ/7772777)
-- Solve for x,y,z
AAns = matrix{{-x^2-y^2+z^2+1, -2*y*z-2*x, 2*x*z-2*y}, {-2*y*z+2*x, -x^2+y^2-z^2+1, -2*x*y-2*z}, { 2*x*z+2*y, -2*x*y+2*z,x^2-y^2-z^2+1}};
eqns = flatten entries((x^2+y^2+z^2+1)*submatrix'(M,{3},{3}) - AAns);
sol = flatten entries matrix first solveSystem(eqns)

xAns = sub(-1,ZZ/7772777)
yAns = sub(-1,ZZ/7772777)
zAns = sub(1,ZZ/7772777)
--xAns = realPart sol_0
--yAns = realPart sol_1
--zAns = realPart sol_2

-*
xAns = sol_0
yAns = sol_1
zAns = sol_2
*-
-- Map the general 6R system to the ring with Kuka parameters
R = ZZ/7772777[x, y, z, tx, ty, tz, S1, S2, S3, C1, C2, C3, D1, D2, D3, R1, R2, R3][ct1, ct2, ct3, st1, st2, st3]; -- we also need to add r's.
DH = {{1,0,1}, {0,1,0}, {675,0,0}, {0,0,0}}; --ENTER SPECIFIC DH PARAMETERS HERE IN THE FOLLOWING ORDER: SIN(ALPHA), COS(ALPHA), D, R
Params = {xAns,yAns,zAns,txAns,tyAns,tzAns}; --x, y, z, tx, ty, tz for specific problem
F = constructEquations();
-- Creating the ring map f
L = ZZ/7772777[ct1, ct2, ct3, st1, st2, st3];
--{ct1=>ct1, ct2=>ct2, ct3=>ct3, ct4=>ct4, ct5=>ct5, ct6=>ct6, st1=>st1, st2=>st2, st3=>st3, st4=>st4, st5=>st5, st6=>st6, S1=>DH#0#0, S2=>DH#0#1, S3=>DH#0#2, S4=>DH#0#3,S5=>DH#0#4, S6=>DH#0#5, C1=>DH#1#0, C2=>DH#1#1, C3=>DH#1#2, C4=>DH#1#3,C5=>DH#1#4, C6=>DH#1#5, D1=>DH#2#0, D2=>DH#2#1, D3=>DH#2#2, D4=>DH#2#3,D5=>DH#2#4,D6=>DH#2#5,R1=>DH#3#0, R2=>DH#3#1, R3=>DH#3#2,R4=>DH#3#3,R5=>DH#3#4,R6=>DH#3#5, x=>Params#0, y=>Params#1, z=>Params#2, tx=>Params#3, ty=>Params#4, tz=>Params#5}
unknowns = {ct1, ct2, ct3, st1, st2, st3};
specification = flatten{unknowns,Params,flatten(DH)}
f = map(L,R,specification);
-- Apply the ring map
func = (poly) -> (f(poly));
lst = apply(F,func)

I=ideal(lst);
dim(I)
degree(I)

