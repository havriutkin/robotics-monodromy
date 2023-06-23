loadPackage "NumericalAlgebraicGeometry";

R = QQ[tx, ty, tz, ax, ay, az][ct1, ct2, ct3, ct4, ct5, ct6, st1, st2, st3, st4, st5, st6];

-- Function to generate random DH parameters
generateRandomDhParameters = dof -> (
    alpha = {};
    r = {};
    d = {};
    theta = {};

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
  T = i -> matrix{{cos(theta#i), -sin(theta#i)*cos(alpha#i),  sin(theta#i)*sin(alpha#i), r#i*cos(theta#i)},
                  {sin(theta#i),  cos(theta#i)*cos(alpha#i), -cos(theta#i)*sin(alpha#i), r#i*sin(theta#i)},
                  {0,             sin(alpha#i),               cos(alpha#i),              d#i},
                  {0,             0,                          0,                         1}};

  -- Compute the forward kinematics transformation matrix
  Tf := product(#d, i -> T(i));

  -- Return the result matrix
  return Tf;
);

-- Construct inverse kinematics equations
constructEquations = (A, alpha, r, d) -> (
    ctheta = {ct1, ct2, ct3, ct4, ct5, ct6};
    stheta = {st1, st2, st3, st4, st5, st6};
    
    -- Construct SE(3) with unknowns
    T = i -> matrix{{ctheta#i, -stheta#i*cos(alpha#i),  stheta#i*sin(alpha#i), r#i*ctheta#i},
                    {stheta#i,  ctheta#i*cos(alpha#i), -ctheta#i*sin(alpha#i), r#i*stheta#i},
                    {0,         sin(alpha#i),           cos(alpha#i),          d#i},
                    {0,         0,                      0,                     1}};
                    
    Tf := product(#d, i -> T(i));

    -- Set correspond elements of A and T to be equal
    myEquations = {};
    for i from 0 to 2 do(
        for j from 0 to 2 do( 
            myEquations = append(myEquations, Tf_(i,j) - A_(i,j));
        );
    );
    
    for i from 0 to #ctheta-1 do(
        myEquations = append(myEquations, (ctheta#i)^2 + (stheta#i)^2 -1);
    );
    
    return myEquations;
);

constructParametriseSE = () -> (
    -- SE(3) along X
    transX = matrix{{1, 0, 0, ax}, 
                    {0, (1-tx^2)/(1+tx^2), (-2*tx)/(1+tx^2), 0}, 
                    {0, (2*tx)/(1+tx^2), (1-tx^2)/(1+tx^2), 0}, 
                    {0, 0, 0, 1}};
    
    -- SE(3) along Y
    transY = matrix{{(1-ty^2)/(1+ty^2), 0, (2*ty)/(1+ty^2), 0}, 
                    {0, 1, 0, ay},
                    {(-2*ty)/(1+ty^2), 0, (1-ty^2)/(1+ty^2), 0}, 
                    {0, 0, 0, 1}};
    
    -- SE(3) along Z
    transZ = matrix{{(1-tz^2)/(1+tz^2), (-2*tz)/(1+tz^2), 0, 0}, 
                    {(2*tz)/(1+tz^2), (1-tz^2)/(1+tz^2), 0, 0},
                    {0, 0, 1, az}, 
                    {0, 0, 0, 1}};
                
    result = transX * transY * transZ;
    return result;
);

testPlanar = () -> (
    -- Define DH parameters randomly
    dof := 2;
    alpha = {0.0, 0.0};
    r = {1.0, 1,0};
    d = {0.0, 0.0};
    theta = {pi/4, pi/4};

    dhParams := {alpha, r, d, theta};

    -- Compute the forward kinematics
    result = forwardKinematics(dhParams#0, dhParams#1, dhParams#2, dhParams#3);

    -- Print the result
    << "Result matrix: " << result << endl << endl;
)

main = () -> (
    --testPlanar();
    
    -- Construct inverse kinematics equations
    --myEquations = constructEquations(result, dhParams#0, dhParams#1, dhParams#2);
    --<< "# of Equations: " << #myEquations << endl;
    --sol = solveSystem myEquations;
    --<< "Solution: " << sol << endl;
);

main()


