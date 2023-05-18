# Finite Difference Formulation

This section outline the formulae used within the code in order to carry out the simulation. Whenever a variable or code snippet refers to the formulation section, this is the part of the code being referred to.

## Matrix A
This is the left hand side coefficient matrix of the system of equations used within the simulation.

\f$ A = \left( \frac{\Delta t D}{2 r_i \Delta r} - \frac{\Delta t D}{2 \Delta r^2} \right) C^{n+1}_{i-1} + \left( 1 + \frac{\Delta t D}{\Delta r^2} \right) C^{n+1}_i - \left( \frac{\Delta t D}{2 r_i \Delta r} + \frac{\Delta t D}{2 \Delta r^2} \right) C^{n+1}_{i+1} \f$

## Matrix B
This is the right hand side coefficient matrix of the system of equations used within the simulation.

\f$ B = \left( \frac{\Delta t D}{2 \Delta r^2} - \frac{\Delta t D}{2 r_i \Delta r} \right) C^{n}_{i-1} + \left( 1 - \frac{\Delta t D}{\Delta r^2} \right) C^{n}_i + \left( \frac{\Delta t D}{2 r_i \Delta r} + \frac{\Delta t D}{2 \Delta r^2} \right) C^{n}_{i+1} \f$
