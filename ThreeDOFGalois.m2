needsPackage "NumericalAlgebraicGeometry";
needsPackage "MonodromySolver";

R = CC[tx, ty, tz][ct1, ct2, ct3, st1, st2, st3];
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

-- Construct inverse kinematics equations
constructEquations = (alpha, r, d) -> (
    ctheta := {ct1, ct2, ct3};
    stheta := {st1, st2, st3};
    
    A := {tx, ty, tz};
    
   
    -- Construct SE(3) with unknowns
    T := i -> matrix{{ctheta#i, -stheta#i*cos(alpha#i),  stheta#i*sin(alpha#i), r#i*ctheta#i},
                    {stheta#i,  ctheta#i*cos(alpha#i), -ctheta#i*sin(alpha#i), r#i*stheta#i},
                    {0,         sin(alpha#i),           cos(alpha#i),          d#i},
                    {0,         0,                      0,                     1}};
    
    Tf := product(#d, i -> T(i));

    -- Set correspond elements of A and T to be equal
    myEquations := {};
   
        for j from 0 to 2 do( 
            myEquations = append(myEquations,Tf_(j,3) - A#j);
        );


    
    
    for i from 0 to #alpha-1 do(
    myEquations = append(myEquations, (ctheta#i)^2 + (stheta#i)^2 -1);
    );
    
    return myEquations;
); 


findTranslationSeed = (alpha, r, d, theta) -> (
    
F := forwardKinematics(alpha, r, d, theta);
ax := F_(0,3);
ay := F_(1,3);
az := F_(2,3);

P1 := point {{ax, ay, az_CC}};
P2 := point {{cos(theta#0), cos(theta#1), cos(theta#2), sin(theta#0), sin(theta#1), sin(theta#2)_CC}};

return {P1, P2};
    
    );

   
end 

restart

load "ThreeDOFGalois.m2"
needsPackage "MonodromySolver";
needsPackage "NumericalAlgebraicGeometry";


DH = generateRandomDhParameters(3);
F = constructEquations(DH#0, DH#1, DH#2);
Seed = findTranslationSeed(DH#0, DH#1, DH#2, DH#3);
S = polySystem(F);
V = first monodromySolve(S, Seed#0, {Seed#1}, Verbose=>true);
netList points V.PartialSols --Extracting some data from the solver
monodromyGroup(V.Graph, FileName => "Planar.gp")
matrix V.BasePoint

evaluate(S,Seed#0,Seed#1) --Extracting data from the system
sub(jacobian S, matrix Seed#1)
help monodromySolve





