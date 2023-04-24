"""
Compute density of a substance that flows around a cylinder (in 2D).
Computations are done using Fenics
"""

from __future__ import print_function
import dolfin
import fenics
import numpy as np
import matplotlib.pyplot as plt
import gmsh

# Define 2D geometry
xmax = 5
ymax = 5
R = 1

LL = dolfin.Point(-xmax,-ymax)
UR = dolfin.Point(xmax,ymax)
domain = mshr.Rectangle(LL,UR)  - mshr.Circle( dolfin.Point(0,0),R)

mesh = mshr.generate_mesh(domain,64)

V = dolfin.FunctionSpace(mesh, 'P', 1)

# Define initial condition
x0 = -3
y0 = -0.5
s0 = 0.01

C0exp = fenics.Expression('1./2./pi/s0/s0*exp(-(pow(x[0]-x0,2.)+pow(x[1]-y0,2.))/2./s0/s0)',
                          pi=np.pi, x0=x0,y0=y0,s0=s0,degree=2)

# Define variational problem
C = dolfin.TrialFunction(V)
v = dolfin.TestFunction(V)

## Drift function for computing advective flow. Perform first
## computations in polar coordinates; then convert to Cartesian
class Drift(dolfin.UserExpression):
    def eval(self,values,x):
        r = np.sqrt(x[0]*x[0] + x[1]*x[1])
        drdt = (1.-1./r/r)*x[0]/r
        dthetadt = -(1.+1./r/r)*x[1]/r/r
        dxdt = drdt*x[0]/r - x[1] * dthetadt
        dydt = drdt*x[1]/r + x[0] * dthetadt 

        values[0] = dxdt
        values[1] = dydt
    def value_shape(self):
        return (2,)


f = Drift()

## New and old solution for time stepping
Csol = dolfin.Function(V)
Cold = dolfin.Function(V)

dt = 0.25

vtkfile = dolfin.File('FwdKolmo-Cylinder/C.pvd')
    

for Pe in [1,10,100]:
    D=1/Pe
    
    # Define initial value
    Cold.assign(fenics.interpolate(C0exp, V))

    ## Define variational problem for updating. This is an implicit
    ## Euler method.
    F = (C-Cold)*v/dt*dolfin.dx + dolfin.dot(D*dolfin.grad(C) - f*C, dolfin.grad(v))*dolfin.dx 

    ## Get matrix formulation of the update
    a, L = dolfin.lhs(F), dolfin.rhs(F)

    num_steps = 10
    
    # Time-stepping
    t = 0
    for n in range(num_steps):
    
        # Update current time
        t += dt
    
        # Compute solution
        dolfin.solve(a == L, Csol)
    
        vtkfile << (Csol,t)
    
        # Update previous solution
        Cold.assign(Csol)
    
    dolfin.plot(Csol,cmap='gray')
    
    plt.savefig("C-"+str(Pe)+".pdf", bbox_inches='tight')

## Create dummy figure for use with makefile to indicate that the program has run
plt.savefig("FwdKolmo-Cylinder.pdf")
