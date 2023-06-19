loadPackage "NumericalAlgebraicGeometry";

R = RR[ct1, ct2, ct3, ct4, ct5, ct6, st1, st2, st3, st4, st5, st6];

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

main = () -> (
    -- Define DH parameters randomly
    dof := 6;
    dhParams := generateRandomDhParameters(dof);

    -- Define the DH parameters for Kuka KR-15/2
    --dof := 6;
    --alpha = {90.0, 0.0, 90.0, -90.0, 90.0, 0.0};
    --r = {0.3, 0.65, 0.155, 0.0, 0.0, 0.0};
    --d = {0.675, 0.0, 0.0, 0.0, 0.0, 0.140};
    --theta = {0.0, -90.0, 0.0, 0.0, 0.0, 0.0};
    --dhParams := {alpha, r, d, theta};

    -- Compute the forward kinematics
    result = forwardKinematics(dhParams#0, dhParams#1, dhParams#2, dhParams#3);

    -- Print the position and orientation
    --<< "DH parameters: " << dhParams << endl << endl;
    --<< "Result matrix: " << result << endl << endl;
    
    -- Construct inverse kinematics equations
    myEquations = constructEquations(result, dhParams#0, dhParams#1, dhParams#2);
    << "# of Equations: " << #myEquations << endl;
    sol = solveSystem myEquations;
    << "Solution: " << sol << endl;
);

main()


