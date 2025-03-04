from dolfin import *
import numpy as np

def init(Name):
    tol = 1E-14 # tolerance for coordinate comparisons  
    
    if Name == 'cantilever':
        Lag,Nx,Ny = [40.0,150,75]
    if Name == 'cantilever_asymmetric':    
        Lag,Nx,Ny = [80.0,150,75]  
    if Name == 'half_wheel':        
        Lag,Nx,Ny = [50.0,200,100]
    if Name == 'bridge':        
        Lag,Nx,Ny = [30.0,200,100]
    if Name == 'MBB_beam':
        Lag,Nx,Ny = [130.0,150,50]                  
    if Name == 'cantilever_twoforces': 
        Lag,Nx,Ny = [60.0,121,121]
    if Name == 'inverter':                            
        Lag,Nx,Ny = [0.01,350,175]
                                   
    if Name == 'cantilever' or Name == 'cantilever_asymmetric':    
        lx,ly = [4.0,1.0]    
        XX,YY = np.meshgrid(np.linspace(0.0,lx,Nx+1),np.linspace(0.0,ly,Ny+1))
        YY -= 0.25 ##shift the pattern by 0.25 units
        #define boundary conditions
        mesh  = RectangleMesh(Point(0.0,0.0),Point(lx,ly),Nx,Ny,'crossed') 
        Vvec = VectorFunctionSpace(mesh, 'CG', 1)         
        class DirBdX(SubDomain):
            def inside(self, x, on_boundary):
                return near(x[0],0.0) 
        class DirBdY(SubDomain):
            def inside(self, x, on_boundary):
                return near(x[0], 0,0) and (0.0 <= x[1] <=0.4)
        dirBdX = DirBdX()
        dirBdY = DirBdY()
        #boundaries = FacetFunction("size_t", mesh)   
        boundaries = MeshFunction("size_t", mesh, mesh.topology().dim()-1, 0)
        boundaries.set_all(0)    
        dirBdX.mark(boundaries, 1)
        dirBdY.mark(boundaries, 2)
        ds = Measure("ds")(subdomain_data=boundaries) 
        zero_x = Constant((0.0, None))
        zero_y = Constant((None, 0.0))
        bcd  = [DircichletBC(Vvec.sub(0), zero_x, boundaries, 1), DirichletBC(Vvec.sub(1),zero_y, boundaries, 2)
            
    if Name == 'cantilever':
        ####################
        # Cantilever case
        ####################
        Load = [Point(lx, 0.5)]
        # Initialize level set function
        #phi_mat = -np.cos(8.0/lx*pi*XX) * np.cos(4.0*pi*YY) - 0.4\
        # + np.maximum(200.0*(0.01-XX**2-(YY-ly/2)**2),.0)\
        # + np.maximum(100.0*(XX+YY-lx-ly+0.1),.0) + np.maximum(100.0*(XX-YY-lx+0.1),.0)    
        phi_mat = -np.cos(3.0/lx*pi*XX) * np.cos(2.0*pi*YY) - 0.6\
         + np.maximum(200.0*(0.01-XX**2-(YY-ly/2)**2),.0)\
         + np.maximum(100.0*(XX+YY-lx-ly+0.1),.0) + np.maximum(100.0*(XX-YY-lx+0.1),.0)                  
                 
    if Name == 'cantilever_asymmetric':    
        ##########################
        # Asymmetric cantilever
        ##########################      
        Load = [Point(lx, 0.75)]   
        # Initialize level set function        
        phi_mat = -np.cos(6.0/lx*pi*XX) * np.cos(4.0*pi*YY) - 0.4\
          + np.maximum(100.0*(XX+YY-lx-ly+0.1),.0)         

    if Name == 'half_wheel' or Name == 'bridge':
        lx,ly = [2.0,1.0]    
        XX,YY = np.meshgrid(np.linspace(0.0,lx,Nx+1),np.linspace(0.0,ly,Ny+1))   
        #define boundary conditions
        mesh  = RectangleMesh(Point(0.0,0.0),Point(lx,ly),Nx,Ny,'crossed') 
        Vvec = VectorFunctionSpace(mesh, 'CG', 1)                
        class DirBd(SubDomain):
            def inside(self, x, on_boundary):
                return abs(x[0])< tol and abs(x[1])< tol 
        class DirBd2(SubDomain):
            def inside(self, x, on_boundary):
                return abs(x[0]-lx)< tol and abs(x[1])< tol
        dirBd,dirBd2 = [DirBd(),DirBd2()]
        boundaries = FacetFunction("size_t", mesh)     
        boundaries.set_all(0)          
        dirBd.mark(boundaries, 1)
        dirBd2.mark(boundaries, 2)  
        ds = Measure("ds")(subdomain_data=boundaries)                     
        Load = [Point(lx/2, 0.0)]

    if Name == 'half_wheel':
        ####################
        # Half-wheel case
        ####################
        bcd  = [DirichletBC(Vvec, (0.0,0.0), dirBd, method='pointwise'),\
               DirichletBC(Vvec.sub(1), 0.0, dirBd2,method='pointwise')]          
        # Initialize level set function
        phi_mat = -np.cos((3.0*pi*(XX-1.0))) * np.cos(7*pi*YY) - 0.3\
         + np.minimum(5.0/ly *(YY-1.0) + 4.0,0) \
         + np.maximum(100.0*(XX+YY-lx-ly+0.1),.0) + np.maximum(100.0*(-XX+YY-ly+0.1),.0)                

    if Name == 'bridge':
        ####################
        # Bridge case
        ####################                    
        bcd  = [DirichletBC(Vvec, (0.0,0.0), dirBd, method='pointwise'),\
               DirichletBC(Vvec, (0.0,0.0), dirBd2,method='pointwise')]                         
        # Initialize level set function
        phi_mat = -np.cos((4.0*pi*(XX-1.0))) * np.cos(4*pi*YY) - 0.2 \
         + np.maximum(100.0*(YY-ly+0.05),.0)  
         
    if Name == 'MBB_beam':
        ####################
        # MBB beam case
        ####################              
        lx,ly = [3.0,1.0]    
        XX,YY = np.meshgrid(np.linspace(0.0,lx,Nx+1),np.linspace(0.0,ly,Ny+1))   
        #define boundary conditions
        mesh  = RectangleMesh(Point(0.0,0.0),Point(lx,ly),Nx,Ny,'crossed') 
        Vvec = VectorFunctionSpace(mesh, 'CG', 1)                
        class DirBd(SubDomain):
            def inside(self, x, on_boundary):
                return near(x[0],.0)   
        class DirBd2(SubDomain):
            def inside(self, x, on_boundary):
                return abs(x[0]-lx) < tol and abs(x[1]) < tol    
        dirBd,dirBd2 = [DirBd(),DirBd2()]
        boundaries = FacetFunction("size_t", mesh)     
        boundaries.set_all(0)          
        dirBd.mark(boundaries, 1)
        dirBd2.mark(boundaries, 2)  
        ds = Measure("ds")(subdomain_data=boundaries)
        bcd  = [DirichletBC(Vvec.sub(0), 0.0, boundaries, 1),\
                DirichletBC(Vvec.sub(1), 0.0, dirBd2,method='pointwise')]                        
        Load = [Point(0.0, 1.0)] 
        # Initialize level set function        
        phi_mat = -np.cos(4.0/lx*pi*XX) * np.cos(4.0*pi*YY) - 0.4\
          + np.maximum(100.0*(XX+YY-lx-ly+0.1),.0)+ np.minimum(5.0/ly *(YY-1.0) + 4.0,0)     
              
    if Name == 'cantilever_twoforces':
        #######################
        # Cantilever two forces
        #######################                         
        lx,ly = [1.0,1.0]    
        XX,YY = np.meshgrid(np.linspace(0.0,lx,Nx+1),np.linspace(0.0,ly,Ny+1))   
        #define boundary conditions
        mesh  = RectangleMesh(Point(0.0,0.0),Point(lx,ly),Nx,Ny,'crossed') 
        Vvec = VectorFunctionSpace(mesh, 'CG', 1)   
        class DirBd(SubDomain):
            def inside(self, x, on_boundary):
                return near(x[0],.0)      
        dirBd = DirBd()
        boundaries = FacetFunction("size_t", mesh)   
        boundaries.set_all(0)    
        dirBd.mark(boundaries, 1)
        ds = Measure("ds")(subdomain_data=boundaries)  
        bcd  = [DirichletBC(Vvec, (0.0,0.0), boundaries, 1)]                         
        Load = [Point(lx, 0.0),Point(lx, 1.0)] 
        # Initialize level set function        
        phi_mat = -np.cos(4.0*pi*(XX-0.5)) * np.cos(4.0*pi*(YY-0.5)) - 0.6 \
        - np.maximum(50.0*(YY-ly+0.1),.0)- np.maximum(50.0*(-YY+0.1),.0)                                                               

    if Name == 'inverter':
        #######################
        # Inverter
        #######################                                          
        lx,ly = [2.0,1.0]    
        XX,YY = np.meshgrid(np.linspace(0.0,lx,Nx+1),np.linspace(0.0,ly,Ny+1))   
        #define boundary conditions
        mesh  = RectangleMesh(Point(0.0,0.0),Point(lx,ly),Nx,Ny,'crossed') 
        Vvec = VectorFunctionSpace(mesh, 'CG', 1)         
        class DirBd(SubDomain):
            def inside(self, x, on_boundary):
                return near(x[0],.0)
                #return (abs(x[0]) < tol and abs(x[1]) < tol) or (abs(x[0]) < tol and abs(x[1]-1.0) < tol) 
        class InputBd(SubDomain):
            def inside(self, x, on_boundary):
                return near(x[0],2.5) and between(x[1], (.47,1))
        class OutputBd(SubDomain):
            def inside(self, x, on_boundary):
                return between(x[0],(3.5,4)) and between(x[1], (0,0.25))
        dirBd = DirBd()
        outputBd = OutputBd()
        inputBd = InputBd()
        
        boundaries = MeshFunction("size_t", mesh, mesh.topology().dim()-1, 0)
        #boundaries = FacetFunction("size_t", mesh)   
        boundaries.set_all(0)    
        dirBd.mark(boundaries, 1)
        outputBd.mark(boundaries, 2)
        inputBd.mark(boundaries, 3)
        ds = Measure("ds")(subdomain_data=boundaries)  
        #bcd  = [DirichletBC(Vvec, (0.0,0.0), dirBd, method='topological')]          
        bcd  = [DirichletBC(Vvec, (0.0,0.0), boundaries, 1)]          
        #Load = [Point(2.5, 0.75)]
        Load = [Point(3.75, 0.15)]
        # Initialize level set function    
        phi_mat = -np.cos(10.0/lx*pi*XX) * np.cos(12.0*pi*YY) - 0.3\
         + np.maximum(100.0*(XX+YY-lx-ly+0.1),.0) + np.maximum(100.0*(XX-YY-lx+0.1),.0)           
             
    return Lag,Nx,Ny,lx,ly,Load,Name,ds,bcd,mesh,phi_mat,Vvec
