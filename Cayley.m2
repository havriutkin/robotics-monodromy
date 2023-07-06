needsPackage "NumericalAlgebraicGeometry";
needsPackage "MonodromySolver";

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

-- Construct inverse kinematics equations
constructEquations = (alpha, r, d) -> (
    ctheta := {ct1, ct2, ct3, ct4, ct5, ct6};
    stheta := {st1, st2, st3, st4, st5, st6};
    
    A := matrix{{x^4+2*x^2*y^2+y^4+2*x^2*z^2+2*y^2*z^2+z^4-6*x^2-6*y^2+2*z^2+1,4*x^3+4*x*y^2+4*x*z^2-8*y*z-4*x, 4*x^2*y+4*y^3+4*y*z^2+8*x*z-4*y, tx*(x^4+2*x^2*y^2+y^4+2*x^2*z^2+2*y^2*z^2+z^4+2*x^2+2*y^2+2*z^2+1)}, {-4*x^3-4*x*y^2-4*x*z^2-8*y*z+4*x,x^4+2*x^2*y^2+y^4+2*x^2*z^2+2*y^2*z^2+z^4-6*x^2+2*y^2-6*z^2+1, 4*x^2*z+4*y^2*z+4*z^3-8*x*y-4*z, ty*(x^4+2*x^2*y^2+y^4+2*x^2*z^2+2*y^2*z^2+z^4+2*x^2+2*y^2+2*z^2+1)}, {-4*x^2*y-4*y^3-4*y*z^2+8*x*z+4*y, -4*x^2*z-4*y^2*z-4*z^3-8*x*y+4*z, x^4+2*x^2*y^2+y^4+2*x^2*z^2+2*y^2*z^2+z^4+2*x^2-6*y^2-6*z^2+1, tz*(x^4+2*x^2*y^2+y^4+2*x^2*z^2+2*y^2*z^2+z^4+2*x^2+2*y^2+2*z^2+1)}, {0,0,0,(x^4+2*x^2*y^2+y^4+2*x^2*z^2+2*y^2*z^2+z^4+2*x^2+2*y^2+2*z^2+1)}};
    
    
    -- Construct SE(3) with unknowns
    T := i -> matrix{{ctheta#i, -stheta#i*cos(alpha#i),  stheta#i*sin(alpha#i), r#i*ctheta#i},
                    {stheta#i,  ctheta#i*cos(alpha#i), -ctheta#i*sin(alpha#i), r#i*stheta#i},
                    {0,         sin(alpha#i),           cos(alpha#i),          d#i},
                    {0,         0,                      0,                     1}};
    
    Tf := product(#d, i -> T(i));

    -- Set correspond elements of A and T to be equal
    myEquations := {};
    for i from 0 to 2 do(
        for j from 0 to 3 do( 
            myEquations = append(myEquations, (x^4+2*x^2*y^2+y^4+2*x^2*z^2+2*y^2*z^2+z^4+2*x^2+2*y^2+2*z^2+1)*Tf_(i,j) - A_(i,j));
        );
    );

    
    
    for i from 0 to #alpha-1 do(
    myEquations = append(myEquations, (ctheta#i)^2 + (stheta#i)^2 -1);
    );
    
    return myEquations;
); 

--Need to find another way to get a seed, now that the system is linearly parametrized 

findSeed = (alpha, r, d, theta) -> (
    
F  := forwardKinematics(alpha, r, d, theta); --Need to multiply this through by denominator term.
A  := matrix{{x^4+2*x^2*y^2+y^4+2*x^2*z^2+2*y^2*z^2+z^4-6*x^2-6*y^2+2*z^2+1,4*x^3+4*x*y^2+4*x*z^2-8*y*z-4*x, 4*x^2*y+4*y^3+4*y*z^2+8*x*z-4*y, tx*(x^4+2*x^2*y^2+y^4+2*x^2*z^2+2*y^2*z^2+z^4+2*x^2+2*y^2+2*z^2+1)}, {-4*x^3-4*x*y^2-4*x*z^2-8*y*z+4*x,x^4+2*x^2*y^2+y^4+2*x^2*z^2+2*y^2*z^2+z^4-6*x^2+2*y^2-6*z^2+1, 4*x^2*z+4*y^2*z+4*z^3-8*x*y-4*z, ty*(x^4+2*x^2*y^2+y^4+2*x^2*z^2+2*y^2*z^2+z^4+2*x^2+2*y^2+2*z^2+1)}, {-4*x^2*y-4*y^3-4*y*z^2+8*x*z+4*y, -4*x^2*z-4*y^2*z-4*z^3-8*x*y+4*z, x^4+2*x^2*y^2+y^4+2*x^2*z^2+2*y^2*z^2+z^4+2*x^2-6*y^2-6*z^2+1, tz*(x^4+2*x^2*y^2+y^4+2*x^2*z^2+2*y^2*z^2+z^4+2*x^2+2*y^2+2*z^2+1)}, {0,0,0,(x^4+2*x^2*y^2+y^4+2*x^2*z^2+2*y^2*z^2+z^4+2*x^2+2*y^2+2*z^2+1)}};
 Equations := {};
    for i from 0 to 2 do(
        for j from 0 to 3 do( 
            Equations = append(Equations, (x^4+2*x^2*y^2+y^4+2*x^2*z^2+2*y^2*z^2+z^4+2*x^2+2*y^2+2*z^2+1)*F_(i,j) - A_(i,j));
        );
    );

s=solveSystem Equations;

return s;

    
    );

   
end 

restart

load "Cayley.m2"
needsPackage "MonodromySolver";
needsPackage "NumericalAlgebraicGeometry";


DH=generateRandomDhParameters(6);
S=constructEquations(DH#0, DH#1, DH#2);
Seed=findSeed(DH#0, DH#1, DH#2, DH#3); --Runs for a long time. Anything we can do about that?



