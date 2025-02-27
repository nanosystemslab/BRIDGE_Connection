# ----------------------------------------------------------------------
# FEniCS 2017.1 code for level set-based structural optimization.
# Written by Antoine Laurain, 2017
# ----------------------------------------------------------------------
from dolfin import *
#parameters["num_threads"] = 1
from init import *
from matplotlib import cm,pyplot as pp
import numpy as np, sys, os
pp.switch_backend('Agg')
#set_log_level(ERROR)
from petsc4py import PETSc
from dolfin import PETScLUSolver

# ----------------------------------------------------------------------
def _main():
    Lag,Nx,Ny,lx,ly,Load,Name,ds,bcd,mesh,phi_mat,Vvec=init(sys.argv[1])      
    eps_er, E, nu = [0.001, 1.0, 0.3]  # Elasticity parameters
    mu,lmbda = Constant(E/(2*(1 + nu))),Constant(E*nu/((1+nu)*(1-2*nu)))
    # Create folder for saving files
    rd = os.path.join(os.path.dirname(__file__),\
     Name +'/LagVol=' +str(np.int_(Lag))+'_Nx='+str(Nx))
    if not os.path.isdir(rd): os.makedirs(rd) 
    # Line search parameters
    beta0_init,ls,ls_max,gamma,gamma2 = [0.5,0,3,0.8,0.8]   
    beta0 = beta0_init
    beta  = beta0
    # Stopping criterion parameters    
    ItMax,It,stop = [int(1.5*Nx), 0, False] 
    # Cost functional and function space
    J = np.zeros( ItMax )  
    V = FunctionSpace(mesh, 'CG', 1)
    VolUnit = project(Expression('1.0',degree=2),V) # to compute volume 
    U = [0]*len(Load) # initialize U
    # Get vertices coordinates 
    gdim     = mesh.geometry().dim()
    dofsV    = V.tabulate_dof_coordinates().reshape((-1, gdim))    
    dofsVvec = Vvec.tabulate_dof_coordinates().reshape((-1, gdim))     
    px,py    = [(dofsV[:,0]/lx)*2*Nx, (dofsV[:,1]/ly)*2*Ny]
    pxvec,pyvec = [(dofsVvec[:,0]/lx)*2*Nx, (dofsVvec[:,1]/ly)*2*Ny]  
    dofsV_max, dofsVvec_max =((Nx+1)*(Ny+1) + Nx*Ny)*np.array([1,2])    
    # Initialize phi  
    phi = Function( V )  
    phi = _comp_lsf(px,py,phi,phi_mat,dofsV_max) 
    # Define Omega = {phi<0}     
    class Omega(SubDomain):
        def inside(self, x, on_boundary):
            return .0 <= x[0] <= lx and .0 <= x[1] <= ly and phi(x) < 0             
    #domains = CellFunction("size_t", mesh)     
    domains= MeshFunction("size_t", mesh, mesh.topology().dim(), 0)
    dX = Measure('dx') 
    n  = FacetNormal(mesh) 
    # Define solver to compute descent direction th
    theta,xi = [TrialFunction(Vvec), TestFunction( Vvec)]     

    # Create an empty PETScMatrix
    petsc_av = PETScMatrix()

    # Assemble the system matrix directly into the PETScMatrix
    assemble((inner(grad(theta), grad(xi)) + 0.1 * inner(theta, xi)) * dX
             + 1.0e4 * (inner(dot(theta, n), dot(xi, n)) * (ds(0) + ds(1) + ds(2))), tensor=petsc_av)

    # Initialize PETScLUSolver with the PETSc matrix
    solverav = PETScLUSolver(petsc_av)

    # Access the underlying PETSc KSP solver
    ksp = solverav.ksp()
    # Access and print information about the preconditioner (PC)
    pc = ksp.getPC()
    # Set KSP and PC options for reuse and performance optimization
    ksp.setType("preonly")
    pc.setType("lu")

    # Setting reuse fill pattern options (if applicable)
    try:
        pc.setReuseFill(True)  # Reuse the fill pattern for efficiency
        pc.setFactorReuseOrdering(True)  # Reuse ordering during factorization
    except AttributeError:
        # In case these settings are not available, catch the error and continue
        print("Factor reuse options are not available in this PETSc version.")

    # Alternatively, set PETSc options globally
    PETSc.Options().setValue("pc_factor_reuse_fill", True)
    PETSc.Options().setValue("pc_factor_reuse_ordering", True)

    # View the current PETSc configuration (for debugging)
    opts = PETSc.Options()
    opts.view()


#---------- MAIN LOOP ----------------------------------------------
    while It < ItMax and stop == False:
        # Update and tag Omega = {phi<0}, then solve elasticity system.  
        omega = Omega()
        domains.set_all(0)
        omega.mark(domains, 1)
        dx = Measure('dx')(subdomain_data = domains)   
        for k in range(0,len(Load)):   
            U[k] = _solve_pde(Vvec,dx,ds,eps_er,bcd,mu,lmbda,Load[k])      
        # Update cost functional 
        compliance = 0
        for u in U:
            eU = sym(grad(u))
            S1 = 2.0*mu*inner(eU,eU) + lmbda*tr(eU)**2
            compliance += assemble( eps_er*S1* dx(0) + S1*dx(1) )  
        vol = assemble( VolUnit*dx(1) )          
        J[It]  = compliance + Lag * vol                            
        # ------- LINE SEARCH ------------------------------------------
        if It > 0 and J[It] > J[It-1] and ls < ls_max:
            ls   += 1
            beta *= gamma            
            phi_mat,phi = [phi_mat_old,phi_old]
            phi_mat = _hj(th_mat, phi_mat, lx,ly,Nx, Ny, beta)
            phi     = _comp_lsf(px,py,phi,phi_mat,dofsV_max) 
            print('Line search iteration : %s' % ls)           
        else:          
            print('************ ITERATION NUMBER %s' % It)                   
            print('Function value        : %.2f' % J[It])
            print('Compliance            : %.2f' % compliance)
            print('Volume fraction       : %.2f' % (vol/(lx*ly))) 
            # Decrease or increase line search step
            if ls == ls_max: beta0 = max(beta0 * gamma2, 0.1*beta0_init)  
            if ls == 0:      beta0 = min(beta0 / gamma2, 1)
            # Reset beta and line search index    
            ls,beta,It = [0,beta0, It+1]
            # Compute the descent direction th           
            th = _shape_der(Vvec,U,eps_er,mu,lmbda,dx,solverav,Lag)
            #th_array = th.vector().array() 
            th_array = th.vector().get_local()
            th_mat = [np.zeros((Ny+1,Nx+1)),np.zeros((Ny+1,Nx+1))]          
            for dof in range(0, dofsVvec_max,2):
                if np.rint(pxvec[dof]) %2 == .0:
                    cx,cy= np.int_(np.rint([pxvec[dof]/2,pyvec[dof]/2]))
                    th_mat[0][cy,cx] = th_array[dof]
                    th_mat[1][cy,cx] = th_array[dof+1]                           
            # Update level set function phi using descent direction th
            phi_old, phi_mat_old = [phi, phi_mat]            
            phi_mat = _hj(th_mat, phi_mat, lx,ly,Nx,Ny, beta)
            if np.mod(It,5) == 0: phi_mat = _reinit(lx,ly,Nx,Ny,phi_mat)     
            phi = _comp_lsf(px,py,phi,phi_mat,dofsV_max)                      
            #------------ STOPPING CRITERION ---------------------------
            if It>20 and max(abs(J[It-5:It]-J[It-1]))<2.0*J[It-1]/Nx**2: 
                stop = True  
            #------------ Plot Geometry --------------------------------  
            if np.mod(It, 10) == 0 or It == 1 or It == ItMax or stop == True:
                # Create a new figure and axis
                fig, ax = pp.subplots()

                # Plot the contourf with appropriate extent and colormap
                contour = ax.contourf(phi_mat, levels=[-10.0, 0.0], extent=[0.0, lx, 0.0, ly], cmap=cm.get_cmap('bone'))

                # Set aspect ratio to equal and box style
                ax.set_aspect('equal', 'box')

                # Add colorbar to make the visualization more informative
                fig.colorbar(contour, ax=ax, orientation='vertical')

                # Set labels and title to make the plot more descriptive
                ax.set_xlabel('X-axis (units)')
                ax.set_ylabel('Y-axis (units)')
                ax.set_title(f'Contour Plot at Iteration {It}')

                # Save the figure before showing it
                pp.savefig(f'output_iteration_{It}.pdf', bbox_inches='tight')

            #if np.mod(It,10)==0 or It==1 or It==ItMax or stop==True:   
            #    # pp.close()     
            #    pp.contourf(phi_mat,[-10.0,.0],extent = [.0,lx,.0,ly],\
            #     cmap=cm.get_cmap('bone'))
            #    pp.axes().set_aspect('equal','box')
            #    # pp.show()
            #    pp.savefig('output.pdf',bbox_inches='tight')             
    return
# ----------------------------------------------------------------------
def _solve_pde(V, dx, ds, eps_er, bcd, mu, lmbda, Load):
    u,v = [TrialFunction(V), TestFunction(V)]
    S1 = 2.0*mu*inner(sym(grad(u)),sym(grad(v))) + lmbda*div(u)*div(v)
    A = assemble( S1*eps_er*dx(0) + S1*dx(1) )   
    b = assemble( inner(Expression(('0.0', '0.0'),degree=2) ,v) * ds(2))    
    U = Function(V)
    delta = PointSource(V.sub(1), Load, -1.0)
    delta.apply(b) 
    for bc in bcd: bc.apply(A,b)    
    # Create an empty PETScMatrix
    #petsc_A = PETScMatrix()
    #solverA = PETScLUSolver(petsc_A(A,b))
    #solver = LUSolver(A)
    #solver.solve(U.vector(), b)    
    # PETSc solver setup
    solver = PETScLUSolver()
    solver.solve(A, U.vector(), b)
    return U
#-----------------------------------------------------------------------    
def _shape_der(Vvec, u_vec , eps_er, mu, lmbda, dx, solver, Lag):   
    xi = TestFunction(Vvec)  
    rv = 0.0 
    for u in u_vec:       
        eu,Du,Dxi = [sym(grad(u)),grad(u),grad(xi)]
        S1 = 2*mu*(2*inner((Du.T)*eu,Dxi) -inner(eu,eu)*div(xi))\
         + lmbda*(2*inner( Du.T, Dxi )*div(u) - div(u)*div(u)*div(xi) )
        rv += -assemble(eps_er*S1*dx(0) + S1*dx(1) + Lag*div(xi)*dx(1))
    th = Function(Vvec)                  
    if isinstance(rv, (float, int)):
        rv_vec = PETScVector()
        rv_vec[:] = rv
        rv = rv_vec
    solver.solve(th.vector(), rv)
    return th
#-----------------------------------------------------------------------        
def _hj(v,psi,lx,ly,Nx,Ny,beta): 
    for k in range(10):
        Dym = Ny*np.repeat(np.diff(psi,axis=0),[2]+[1]*(Ny-1),axis=0)/ly 
        Dyp = Ny*np.repeat(np.diff(psi,axis=0),[1]*(Ny-1)+[2],axis=0)/ly
        Dxm = Nx*np.repeat(np.diff(psi),[2]+[1]*(Nx-1),axis=1)/lx 
        Dxp = Nx*np.repeat(np.diff(psi),[1]*(Nx-1)+[2],axis=1)/lx          
        g = 0.5*( v[0]*(Dxp + Dxm) + v[1]*(Dyp + Dym)) \
          - 0.5*(np.abs(v[0])*(Dxp - Dxm) + np.abs(v[1])*(Dyp - Dym)) 
        maxv = np.max(abs(v[0]) + abs(v[1]))
        dt  = beta*lx / (Nx*maxv)
        psi = psi - dt*g
    return  psi         
#-----------------------------------------------------------------------         
def _reinit(lx,ly,Nx,Ny,psi):      
    Dxs = Nx*(np.repeat(np.diff(psi),[2]+[1]*(Nx-1),axis=1) \
          +np.repeat(np.diff(psi),[1]*(Nx-1)+[2],axis=1))/(2*lx) 
    Dys = Ny*(np.repeat(np.diff(psi,axis=0),[2]+[1]*(Ny-1),axis=0)\
          +np.repeat(np.diff(psi,axis=0),[1]*(Ny-1)+[2],axis=0))/(2*ly)                  
    signum = psi / np.power(psi**2 + ((lx/Nx)**2)*(Dxs**2+Dys**2),0.5)
    for k in range(0,2):      
        Dym = Ny*np.repeat(np.diff(psi,axis=0),[2]+[1]*(Ny-1),axis=0)/ly 
        Dyp = Ny*np.repeat(np.diff(psi,axis=0),[1]*(Ny-1)+[2],axis=0)/ly
        Dxm = Nx*np.repeat(np.diff(psi),[2]+[1]*(Nx-1),axis=1)/lx 
        Dxp = Nx*np.repeat(np.diff(psi),[1]*(Nx-1)+[2],axis=1)/lx               
        Kp = np.sqrt((np.maximum(Dxm,0))**2 + (np.minimum(Dxp,0))**2 \
           + (np.maximum(Dym,0))**2 + (np.minimum(Dyp,0))**2)
        Km = np.sqrt((np.minimum(Dxm,0))**2 + (np.maximum(Dxp,0))**2 \
           + (np.minimum(Dym,0))**2 + (np.maximum(Dyp,0))**2)           
        g  = np.maximum(signum,0)*Kp + np.minimum(signum,0)*Km
        psi  = psi - (0.5*lx/Nx)*(g - signum)       
    return psi   
#-----------------------------------------------------------------------
def _comp_lsf(px,py,phi,phi_mat,dofsV_max):               
    for dof in range(0,dofsV_max):              
        if np.rint(px[dof]) %2 == .0: 
            cx,cy = np.int_(np.rint([px[dof]/2,py[dof]/2]))                                            
            phi.vector()[dof] = phi_mat[cy,cx]
        else:
            cx,cy = np.int_(np.floor([px[dof]/2,py[dof]/2]))                      
            phi.vector()[dof] = 0.25*(phi_mat[cy,cx] + phi_mat[cy+1,cx]\
              + phi_mat[cy,cx+1] + phi_mat[cy+1,cx+1])    
    return phi            
# ----------------------------------------------------------------------
if __name__ == '__main__':
    _main()

