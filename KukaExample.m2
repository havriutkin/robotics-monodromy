restart

load "6RFinalSymbolic.m2"
load "forwardKinematics.m2"
needsPackage "MonodromySolver"
needsPackage "NumericalAlgebraicGeometry"

S = CC[x,y,z];
--DH parameters for Kuka
alpha = {-pi/2,0,pi/2,-pi/2,pi,0};
d = {675,100,0,600,0,140};
r = {300,650,155,0,0,0};
theta = {0,-pi/2,0,0,0,0};
-- SE(3) matrix for Kuka
M = forwardKinematics(alpha,r,d,theta)
-- x,y,z,tx,ty,tz for Kuka
txAns = sub(-300,ZZ/7772777)
tyAns = sub(-40,ZZ/7772777)
tzAns = sub(1480,ZZ/7772777)
-- Solve for x,y,z
AAns = matrix{{-x^2-y^2+z^2+1, -2*y*z-2*x, 2*x*z-2*y}, {-2*y*z+2*x, -x^2+y^2-z^2+1, -2*x*y-2*z}, { 2*x*z+2*y, -2*x*y+2*z,x^2-y^2-z^2+1}};
eqns = flatten entries((x^2+y^2+z^2+1)*submatrix'(M,{3},{3}) - AAns);
sol = flatten entries matrix first solveSystem(eqns)

xAns = sub(1,ZZ/7772777)
yAns = sub(1,ZZ/7772777)
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
R = ZZ/7772777[x, y, z, tx, ty, tz, S1, S2, S3, S4, S5, S6, C1, C2, C3, C4, C5, C6, D1, D2, D3, D4, D5, D6, R1, R2, R3, R4, R5, R6][ct1, ct2, ct3, ct4, ct5, ct6, st1, st2, st3, st4, st5, st6]; -- we also need to add r's.
DH = {{1,0,1,-1,1,0}, {0,1,0,0,0,1}, {675,0,0,0,600,140}, {300,650,155,0,0,0}}; --ENTER SPECIFIC DH PARAMETERS HERE IN THE FOLLOWING ORDER: SIN(ALPHA), COS(ALPHA), D, R
Params = {xAns,yAns,zAns,txAns,tyAns,tzAns}; --x, y, z, tx, ty, tz for specific problem
F = constructEquations();
-- Creating the ring map f
L = ZZ/7772777[ct1, ct2, ct3, ct4, ct5, ct6, st1, st2, st3, st4, st5, st6];
unknowns = {ct1, ct2, ct3, ct4, ct5, ct6, st1, st2, st3, st4, st5, st6};
specification = flatten{unknowns,Params,flatten(DH)}
f = map(L,R,specification);
-- Apply the ring map
func = (poly) -> (f(poly));
lst = apply(F,func)

I=ideal(lst);
dim(I)
degree(I)

J = ideal(lst_{0,1,2,3,4,5} | take(lst,-6))
dim(J)
degree J
