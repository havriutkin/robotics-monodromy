-*

    This package contains functions related to a standard kinematics of a robot.

    Functions in this package:
        - generateRandomDhParameters (dof: int):
            generates random set of dh parameters with give dof and following boundries:
                alpha from (-pi/2, pi/2)
                r from (-10, 10)
                d from (0, 10)
                theta from (-pi, pi)
            returns list in format [alpha, r, d, theta].
        
        - forwardKinematics (alpha, r, d, theta):
            performs forward kinematics with given DH parameters.
            returns SE(3) matrix which represents EE position.

        - constructEquations (A: SE(3), alpha, r, d, thetea)
            generates and returns list of inverse kinematics equations.
                A - SE(3) which represents EE position (could be parametrised).

*-

-- Function generates random DH parameters
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

-- Funtion performs forward kinematics
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

-- Function constructs inverse kinematics equations
constructEquations = (A, alpha, r, d) -> (
    -- Define local rings
    R := CC[ct1, ct2, ct3, ct4, ct5, ct6, st1, st2, st3, st4, st5, st6];

    ctheta := {ct1, ct2, ct3, ct4, ct5, ct6};
    stheta := {st1, st2, st3, st4, st5, st6};
    
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
            myEquations = append(myEquations, Tf_(i,j) - A_(i,j));
        );
    );
    
    -- Add equations that constraint cos and sin
    for i from 0 to #alpha-1 do(
        myEquations = append(myEquations, (ctheta#i)^2 + (stheta#i)^2 -1);
    );
    
    return myEquations;
);