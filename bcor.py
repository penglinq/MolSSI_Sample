# Define symmetry of correlated bath basis C_EO_LO parameterized by a real 
# zero-diagonal anti-symmetric matrix. Provide methods to parameterize a 
# basis C_lo_eo and to convert parameters to the basis C_lo_eo.
# This class can be auto-differentiated with Jax.
# Reference: Hunter & Parrinello JCP 1994
# Author: Linqing Peng <linqingp@outlook.com>

import itertools as it
import numpy as np
import numpy.linalg as la
import libdmet_solid.utils.logger as log
from libdmet_solid.routine.numpy_helper import expm, block
import jax.numpy as diffnp

class Bcor(object):
    def __init__(self, lattice=None, restricted=True, nbath=None, param=None):
        self.param = param
        if param != None:
            self.value = self.evaluate()
        else:
            self.value = None
        if nbath is None:
            # Assume nbath = # of orbitals in the first cell
            self.nbath = self.nscsites
        else:
            self.nbath = nbath
        if lattice != None:
            self.nenv = lattice.nscsites * (lattice.ncells - 1)
            self.nscsites = lattice.nscsites
        else: 
            self.nenv = None
            self.nscsites = None
        if restricted:
            self.spin = 1
        else: 
            self.spin = 2
    
    def init_from_basis(self, C_lo_eo, nbath=None, lattice=None):  
        '''
            Initialize a bcor object with bcor.param from parametrizing
            an initial guess of the rotation matrix. 

            Args:
                C_lo_eo: the rotation matrix from the localized orbitals (lo)  
                    to the embedding orbitals (eo). 
                    Shape: spin * nscsites * nlo * neo
                nbath: the number of bath orbitals
                lattice: the k-symmetric psycf lattice calculated 
        '''
        log.eassert(C_lo_eo.ndim == 4, "dimension of input basis is not 4")
        if nbath is None:
            # Assume nbath = # of orbitals in the first cell
            self.nbath = self.nscsites
        else:
            self.nbath = nbath
        if lattice != None:
            self.nenv = lattice.nscsites * (lattice.ncells - 1)
            self.nscsites = lattice.nscsites
        self.spin = C_lo_eo.shape[0]
        x = np.zeros((self.spin, self.nbath, self.nenv - self.nbath))                  
        for s in range(self.spin):
            C = C_lo_eo[s,1:,:,self.nscsites:].reshape(self.nenv, self.nbath)
            C1 = C[:self.nbath]
            C2 = C[self.nbath:]
            v, sigma, wt = la.svd(C1) 
            U = np.dot(wt.T, v.T)
            C1_tilde = np.dot(C1, U)
            C2_tilde = np.dot(C2, U)
            v1, sigma1, wt1 = la.svd(C1_tilde)
            P_sqrt = np.dot(v1, np.dot(np.diag(np.arccos(sigma1)), wt1))
            # Note: Corrected a mistake in Hunter1994 in the expression of X.
            sin_P_sqrt = np.dot(v1, np.dot(np.diag(np.sin(np.arccos(sigma1))), wt1))
            x[s] = -np.dot(P_sqrt, np.dot(la.inv(sin_P_sqrt), C2_tilde.T))
        self.param = np.ndarray.flatten(x)
        self.update(self.param)

        # Check symmetry
        # Check if the parameterized (symmetrized) basis matches with the 
        # input C_lo_eo
        Cin = C_lo_eo.reshape(self.spin, self.nenv + self.nscsites, \
                self.nscsites + self.nbath) 
        Ccalc = self.get().reshape(self.spin, self.nenv + self.nscsites, \
                self.nscsites + self.nbath)
        Ddiff = sum([la.norm(Cin[s].dot(Cin[s].T) - Ccalc[s].dot(Ccalc[s].T)) \
                for s in range(self.spin)])
        log.check(Ddiff < 1e-7, "symmetrization imposed on initial guess, \
                    is not satisfied")
        
    def update(self, param, autodiff=False):
        '''
            Update the paramters of a bcor object. 

            Arg:
                param: The new parameters. 1d numpy array.
                autodiff: if True, use the jax.numpy to carry out 
                          auto-differentiation. Otherwise, use the normal numpy
        '''
        if autodiff:
            numpy = diffnp
        else:
            numpy = np
        self.param = numpy.array(param, copy=False)
        self.value = self.evaluate(autodiff=autodiff)

    def get(self, i = 0, kspace = True):
        '''
            Return C_lo_eo to which Bcor corresponds.

            Return:
                The rotation matrix from the localized orbitals (lo) to 
                the embedding orbitals (eo). Shape: spin * nscsites * nlo * neo
        '''
        log.eassert(self.value is not None, "Bcor not initialized yet")
        return self.value

    def evaluate(self, autodiff=False):
        '''
            Convert a bcor object to the corresponding rotation matrix from 
            localized orbitals (lo) to the embedding orbitals (eo) C_lo_eo. 
            All the notations follow Hunter & Parrinello JCP 1994.
            
            Return:
                C_lo_eo: the rotation matrix from the localized orbitals (lo)  
                    to the embedding orbitals (eo). 
                    Shape: spin * nscsites * nlo * neo
        '''
        if autodiff:
            numpy = diffnp
        else:
            numpy = np
        ncells = self.nenv // self.nscsites + 1
        x = self.param.reshape((self.spin, self.nbath, self.nenv - self.nbath))
        for s in range(self.spin):
            A = block(numpy.zeros((self.nbath,self.nbath)), x[s], \
                    -x[s].T, numpy.zeros((self.nenv - self.nbath,self.nenv - \
                            self.nbath)), autodiff=autodiff)
            C = expm(A, autodiff=autodiff)[:, :self.nbath]
            C = block(numpy.eye(self.nscsites), numpy.zeros((self.nscsites, \
                            self.nbath)), numpy.zeros((self.nenv, \
                            self.nscsites)), C, autodiff=autodiff)
            if s==0:
                C_lo_eo = C.reshape((ncells, self.nscsites, self.nscsites + \
                            self.nbath))[numpy.newaxis]
            else:
                C_lo_eo = numpy.append(C_lo_eo, C.reshape((ncells, 
                            self.nscsites, self.nscsites + \
                            self.nbath))[numpy.newaxis], axis=0)
        return C_lo_eo

    def assign(self, C_lo_eo):
        '''
            Given a rotation matrix , calculate the corresponding 
            parameters. Update the bcor object with the new parameters.
            All the notations follow Hunter & Parrinello JCP 1994.
            
            Args:
                C_lo_eo: the rotation matrix from the localized orbitals (lo) 
                    to the embedding orbitals (eo). 
                    Shape: spin * nscsites * nlo * neo
        '''
        basis_shape = tuple((self.spin, self.nenv // self.nscsites + 1, \
                    self.nscsites, self.nscsites + self.nbath))
        log.eassert(C_lo_eo.shape == basis_shape, "The correlation potential \
                should have shape %s, rather than %s", \
                basis_shape, C_lo_eo.shape)
        
        x = np.zeros((self.spin, self.nbath, self.nenv - self.nbath))                  
        for s in range(self.spin):
            C = C_lo_eo[s, 1:, :, self.nscsites:].reshape(self.nenv, self.nbath)
            C1 = C[:self.nbath]
            C2 = C[self.nbath:]
            v, sigma, wt = la.svd(C1) 
            U = np.dot(wt.T, v.T)
            C1_tilde = np.dot(C1, U)
            C2_tilde = np.dot(C2, U)
            v1, sigma1, wt1 = la.svd(C1_tilde)
            P_sqrt = np.dot(v1, np.dot(np.diag(np.arccos(sigma1)), wt1))
            sin_P_sqrt = np.dot(v1, np.dot(np.diag(np.sin(np.arccos(sigma1))), wt1))
            x[s] = - np.dot(P_sqrt, np.dot(la.inv(sin_P_sqrt), C2_tilde.T))
        self.param = np.ndarray.flatten(x)
        self.update(self.param)
        
        # Check symmetry
        # check if the parameterized (symmetrized) basis matches with the 
        # input C_lo_eo
        Cin = C_lo_eo.reshape(self.spin, self.nenv + self.nscsites, \
                self.nscsites + self.nbath) 
        Cout = self.get().reshape(self.spin, self.nenv + self.nscsites, \
                self.nscsites + self.nbath)
        Ddiff = sum([la.norm(Cin[s].dot(Cin[s].T) - Cout[s].dot(Cout[s].T)) \
                for s in range(self.spin)])
        log.check(Ddiff < 1e-7, "symmetrization imposed on initial guess, \
                    is not satisfied")
        
    def __str__(self):
        return self.evaluate().__str__()

def BcorLocal(lattice, restricted=True, nbath=None, param=None):
    bcor = Bcor(lattice, restricted, nbath, param)
    return bcor

def test():
    from libdmet_solid.system import lattice
    from pyscf.pbc import gto
    R = 1.00
    cell = gto.Cell()
    cell.a = ''' 10.0    0.0     0.0
                 0.0     10.0    0.0
                 0.0     0.0     %s '''%(R*2.0)
    cell.atom = ''' H 5.0      5.0      0
                    H 5.0      5.0      %s '''%R
    cell.basis = 'gth-dzv'
    cell.precision = 1e-12
    cell.build(unit='Angstrom')
    kmesh = [1, 1, 3]
    Lat = lattice.Lattice(cell, kmesh)
    np.set_printoptions(4, linewidth=1000, suppress=False) 
    
    log.result("Test converting from parameters to C_lo_eo")
    bcor = BcorLocal(Lat, True, 2)
    param0 = np.asarray([-0.0041, -0.002, -0.1222, -0.2423, 0.116, 0.114, \
            0.0711, 0.0355, 0.1932, 0.351, -0.0866, -0.0511])
    bcor.update(param0)
    log.result("C_lo_eo:\n%s", bcor.get())
    
    log.result("Test converting from C_lo_eo to parameters")
    bcor.assign(bcor.get())
    assert np.max(np.abs(param0 - bcor.param)) < 1E-07

if __name__ == "__main__":
    test()
