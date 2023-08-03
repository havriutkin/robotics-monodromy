# Using Monodromy to generate Galois group in 6-R robotic's arm 

This repository contains the code and documentation for an ongoing research project focused on using monodromy to solve inverse kinematic problem and to generate Galois group associated with given robot.

## Abstract

This project aims to develop a novel solution for the inverse kinematics problem in 6R wrist-partitioned robots using monodromy group analysis. The Denavit-Hartenberg parameters are employed to describe the robot's kinematics, and the positions and orientations of the end effector are represented in the SE(3) space. By setting up a parametrized polynomial system and leveraging homotopy continuation techniques, we investigate the solution variety and generate the monodromy group. This group analysis is expected to provide valuable information about the robot's behavior and geometric characteristics, ultimately enhancing its performance in practical applications.

## Current result
At this point, we have implemented everything related to the robot's kinematics. Functions for using monodromy and generating the monodromy group have been created. The next steps are to analyze the generated group.
Work in progress! Bugs can be encountered!
