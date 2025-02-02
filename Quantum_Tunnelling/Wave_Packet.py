from scipy import sparse as sparse
import numpy as np 
from scipy.sparse import linalg as ln


class wave_packet():
    def __init__(self, sigma0 : float = 5.0, x0 : int = -150, k0 : float = 1.0 , x_end : int = 200, x_begin : int = -200,
                  dt : float = 0.01, barrier_height : float = 1.0, barrier_width : float = 4.0, n_points : int = 100):

        '''
        Here in the Constructor, we:
        - Discretize the x-axis
        - Create the Gaussian Wave-packet array. 
        - Create the potential barrier
        - Create the hamiltonian array, using the Finite Difference Approach, and store the result in a sparse matrix
        - We also define the evolution matrix, which we derived in the attached PDF. This is stored in a compressed sparse row format.
        '''
        self.sigma = sigma0
        self.k = k0
        self.dt = dt
        self.barrierheight = barrier_height
        self.barrierwidth = barrier_width
        self.n = n_points
        self.t = 0

        #discretize the x-axis
        self.xgrid, self.dx  = np.linspace(x_begin,x_end,n_points,retstep=True)
        
        # Initialize wave function as a complex array
        norm = (2.0 * np.pi * sigma0**2)**(-0.25)
        self.psi = np.exp(-(self.xgrid - x0)**2 / (4.0 * sigma0**2)).astype(np.complex128)
        self.psi *= np.exp(1.0j * k0 * self.xgrid)
        self.psi *= norm

        #Create the Potential Barrier
        self.potential = np.array([self.barrierheight if 0 <  x < self.barrierwidth else 0 for x in self.xgrid])

        #Create the Hamiltonian using the Finite Difference Method and Sparse Matrices
        h_diag = np.ones(n_points) / self.dx**2 + self.potential
        h_non_diag = np.ones(n_points-1)*(-0.5/self.dx**2)
        self.hamiltonian = sparse.diags([h_non_diag,h_diag,h_non_diag],[-1,0,1])

        #Create the Crank-Nicholsen Time Evolution Step using Implicit + Explicit Euler
        imp = (sparse.eye(n_points) - self.hamiltonian/2*1.0j).astype(np.complex128).tocsc()
        exp = (sparse.eye(n_points) + self.hamiltonian/2*1.0j).astype(np.complex128).tocsc()
        self.evolution_matrix = ln.inv(imp).dot(exp).tocsr()
    
    
    def evolve(self):
        '''
        Using the Crank-Nicolsen Time-Evolution operator to 'evolve' the wavefunction
        '''
        self.psi = self.evolution_matrix.dot(self.psi)
        #The probability density of the new wavefunction
        self.prob = abs(self.psi)**2
        #Normalising both the updated wavefunction and the probability density
        self.norm = sum(self.prob)*self.dx #Discrete Integration so dx neccessary
        self.prob /= self.norm
        self.psi /= self.norm**0.5
        self.t += self.dt 
        return self.prob

    
    def calculate_probabilities(self):
        '''
        Calculating the probability of finding the particle at either side of the barrier, at a time t
        '''
        left_indices = np.where(self.xgrid < 0)
        right_indices = np.where(self.xgrid > self.barrierwidth)

        #Discrete probabilities
        P_left = np.sum(self.prob[left_indices]) * self.dx
        P_right = np.sum(self.prob[right_indices]) * self.dx
        time = self.t
        
        return P_left, P_right, time

    def analytical_transmission_coefficient(self, E : float):
        """
        Compute the analytical transmission coefficient T using the standard quantum mechanical formula. This is derived in the accompanying PDF
        """
        V0 = self.barrierheight
        a = self.barrierwidth
        m = 1  # Assume m = 1 in atomic units where ħ = 1
    
        if E >= V0:
            return 1.0  # Classical behavior: full transmission
    
        k2 = np.sqrt(2 * m * (V0 - E))  # Decay inside barrier
        k1 = np.sqrt(2 * m * E)  # Wavenumber outside barrier
    
        T = 1 / (1 + (V0**2 * np.sinh(k2 * a)**2) / (4 * E * (V0 - E)))
        return T
