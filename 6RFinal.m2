loadPackage "NumericalAlgebraicGeometry";
loadPackage "MonodromySolver";

R = CC[x, y, z, tx, ty, tz][ct1, ct2, ct3, ct4, ct5, ct6, st1, st2, st3, st4, st5, st6];
-- Replaced ci si with the angles themselves. Parametrization is linear, but still we have problems.
-- Function to generate random DH parameters
generateRandomDhParameters = dof -> (
    alpha := {};
    r := {};
    d := {};
    theta := {};

    -- Generate random values
    for i from 1 to dof do (
        alpha = append(alpha, random(-pi/2, pi/2));
        r = append(r, random(-10, 10));
        d = append(d, random(0, 10));
        theta = append(theta, random(-pi, +pi));
    );

    return {alpha, r, d, theta};
);

-- Forward kinematics function
forwardKinematics = (alpha, r, d, theta) -> (
  -- Define the transformation matrix for one joint
  T := i -> matrix{{cos(theta#i), -sin(theta#i)*cos(alpha#i),  sin(theta#i)*sin(alpha#i), r#i*cos(theta#i)},
                   {sin(theta#i),  cos(theta#i)*cos(alpha#i), -cos(theta#i)*sin(alpha#i), r#i*sin(theta#i)},
                   {0,             sin(alpha#i),               cos(alpha#i),              d#i},
                   {0,             0,                          0,                         1}};

  -- Compute the forward kinematics transformation matrix
  Tf := product(#d, i -> T(i));

  -- Return the result matrix
  return Tf;
);



constructEquations = (alpha, r, d) -> (
    ctheta := {ct1, ct2, ct3, ct4, ct5, ct6};
    stheta := {st1, st2, st3, st4, st5, st6};
    
    A := matrix{{-x^2-y^2+z^2+1, -2*y*z-2*x, 2*x*z-2*y}, {-2*y*z+2*x, -x^2+y^2-z^2+1, -2*x*y-2*z}, { 2*x*z+2*y, -2*x*y+2*z,x^2-y^2-z^2+1}};
    B := {tx, ty, tz};
    
    
    
    
    -- Construct SE(3) with unknowns
    T := i -> matrix{{ctheta#i, -stheta#i*cos(alpha#i),  stheta#i*sin(alpha#i), r#i*ctheta#i},
                    {stheta#i,  ctheta#i*cos(alpha#i), -ctheta#i*sin(alpha#i), r#i*stheta#i},
                    {0,         sin(alpha#i),           cos(alpha#i),          d#i},
                    {0,         0,                      0,                     1}};
    
    Tf := product(#d, i -> T(i));

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

    for l from 0 to #alpha-1 do(
    myEquations = append(myEquations, (ctheta#l)^2 + (stheta#l)^2 -1);
    );
    
    return myEquations;

);

findSeed = (alpha, r, d, theta) -> (
    
F := forwardKinematics(alpha, r, d, theta);
ax := F_(0,3);
ay := F_(1,3);
az := F_(2,3);

R := CC[x,y,z];
--Challenging to find seed in general, may want to take wrist partitioned approach! 

S := submatrix'(F, {3}, {3});

A := matrix{{-x^2-y^2+z^2+1, -2*y*z-2*x, 2*x*z-2*y}, {-2*y*z+2*x, -x^2+y^2-z^2+1, -2*x*y-2*z}, { 2*x*z+2*y, -2*x*y+2*z,x^2-y^2-z^2+1}};
  
    Equations := {};
    for i from 0 to 2 do(
        for j from 0 to 2 do( 
            Equations = append(Equations, (x^2+y^2+z^2+1)*S_(i,j) - A_(i,j));
        );
    );

Sols := solveSystem Equations;

R = CC[x, y, z, tx, ty, tz][ct1, ct2, ct3, ct4, ct5, ct6, st1, st2, st3, st4, st5, st6];

P1 := point {{(coordinates(Sols#0))#0, (coordinates(Sols#0))#1, ((coordinates(Sols#0))#2), ax, ay, az_CC}};
P2 := point {{cos(theta#0), cos(theta#1), cos(theta#2),cos(theta#3), cos(theta#4), cos(theta#5), sin(theta#0), sin(theta#1), sin(theta#2),sin(theta#3), sin(theta#4), sin(theta#5)_CC}};

return {P1, P2};
    
    );


end 

restart

load "6RFinal.m2"
needsPackage "MonodromySolver"
needsPackage "NumericalAlgebraicGeometry"
errorDepth=0;

DH=generateRandomDhParameters(6);
G = constructEquations(DH#0, DH#1, DH#2);
Seed = findSeed(DH#0, DH#1, DH#2, DH#3);
--Sys = polySystem(G);
Sys = polySystem(G_{0,1,2,3,4,5,6,7,8,9,10,11}|take G_(0,-6)); --Get full information from "relaxation"! Need to find a way to specialize system in nonsingular way, and we get full Galois group info.
--evaluate(Sys, Seed#0, Seed#1)
--The following code finds a suitable "relaxation" of the system, by cycling through (12) choose (6) different ways (6-subsets to remove) to remove equations and then taking the product of the evaluation of SVD at jacobian and comparing to zero to see if it is zero or not. 




first SVD(sub(jacobian Sys, matrix Seed#1)) --If product of SVD decomp is near 0, likely singular
V = first monodromySolve(Sys, Seed#0, {Seed#1}, Verbose=>true);  
netList points V.PartialSols --Extracting some data from the solver
monodromyGroup(V.Graph, FileName => "ThreeWrist.gp")
matrix V.BasePoint

--Find way to compute best "relaxed" systems
--Beautify code in symbolic computation to compute degree 







    
