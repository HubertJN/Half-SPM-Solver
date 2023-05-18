# Finite Difference Formulation

This section outline the formulae used within the code in order to carry out the simulation. Whenever a variable or code snippet refers to the formulation section, this is the part of the code being referred to.

## Matrix A
This is the left hand side coefficient matrix of the system of equations used within the simulation, of size \f$ n \times n \f$ with elements given by

\f[
A_{i,j} = \begin{cases}
          - \frac{\Delta t D}{2 \Delta r^2} + \frac{\Delta t D}{2 r_i \Delta r} & j = i-1 \text{ for } 2 \leqslant i \leqslant n-1\\ 
          1 + \frac{\Delta t D}{\Delta r^2} & j = i \\ 
          - \frac{\Delta t D}{2 r_i \Delta r} - \frac{\Delta t D}{2 \Delta r^2} & j = i+1 \text{ for } 2 \leqslant i \leqslant n-1\\ 
          - \frac{\Delta t D}{\Delta r^2} & \left(i,j\right) = \left(1,2\right), \left(n,n-1\right) 
          \end{cases} 
\f]

## Matrix B
This is the right hand side coefficient matrix of the system of equations used within the simulation, of size \f$ n \times n \f$ with elements given by

\f[
B_{i,j} = \begin{cases}
          \frac{\Delta t D}{2 \Delta r^2} - \frac{\Delta t D}{2 r_i \Delta r} & j = i-1 \text{ for } 2 \leqslant i \leqslant n-1\\ 
          1 - \frac{\Delta t D}{\Delta r^2} & j = i \\ 
          \frac{\Delta t D}{2 r_i \Delta r} + \frac{\Delta t D}{2 \Delta r^2} & j = i+1 \text{ for } 2 \leqslant i \leqslant n-1\\ 
          \frac{\Delta t D}{\Delta r^2} & \left(i,j\right) = \left(1,2\right), \left(n,n-1\right) 
          \end{cases} 
\f]
