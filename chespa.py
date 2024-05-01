import numpy as np
'''
equations are from this paper and supplementary:

Boulton, S., Akimoto, M., Selvaratnam, R., Bashiri, A., & Melacini, G. (2014). 
A tool set to map allosteric networks through the NMR chemical shift covariance analysis. 
Scientific Reports, 4, 7306.
'''
class Chespa:
    '''
    Chemical shift projection analysis
    Useful in revealing false positive clustering
    
    Takes the HSQC 2D chemical shift values from a reference state (ref),
    agonist/activated state (B), and 
    an antagonist/reverse agonist state (A)
    '''
    def __init__(self,
                 file,
                 het_nuc='N'):
        '''
        The input csv file must have columns labeled {state}w1, {state}w2 
        for states ref, A, and B
        '''
        df = pd.read_csv(file)
        if het_nuc == 'N':
            het_coef = 0.2
        if het_nuc == 'C':
            het_coef = 0.25
        # hetero nucleus shifts in column 0, H shifts in column 1 
        self.states = {
        'ref': df[['refw1','refw2']].values,
        'A': df[['Aw1','Aw2']].values,
        'B': df[['Bw1','Bw2']].values
        }
        self.ref_to_B = self.get_state_distance('B')
        self.ref_to_A = self.get_state_distance('A')
        self.state_angles = self.get_state_angle()

    def get_state_distance(self, state):
        '''
        SI eq 5
        magnitude of vector between ref and state
        '''
        
        Hrb= (ref[:,1]-self.states[state][:,1])**2
        Hetrb = (het_coef*(ref[:,0]-self.states[state][:,0]))**2
        return np.sqrt(Hrb+Hetrb)
        
    def get_state_angle(self):
        '''
        SI eq 6
        '''
        ref, A, B = list(self.states.values())
        # normalized vectors from ref
        vb = (B-ref)/np.linalg.norm((B-ref))
        va = (A-ref)/np.linalg.norm((A-ref))
        # row wise dot product of 2 n_residues X 2 arrays
        dot = (va*vb).sum(axis=1)
        # multiply magnitudes
        denominator = (self.ref_to_A*self.ref_to_B) 
        cos_theta = dot/denominator
        # fix divide by 0's
        cos_theta = np.nan_to_num(cos_theta, nan=0)
        return cos_theta
    
    def get_fractional_activation(self):
        '''
        SI eq 7
        Normalized projection of antagonist A onto agonist B
        from text:
        The fractional activation reveals whether a given perturbation 
        shifts the protein  towards activation (X > 0) or inactivation (X < 0).
        '''
        ref, A, B = list(self.states.values())
        # attribute?
        # normalized vectors from ref
        vb = (B-ref)/np.linalg.norm((B-ref))
        va = (A-ref)/np.linalg.norm((A-ref))
        # row wise dot product of 2 n_residues X 2 arrays
        dot = (va*vb).sum(axis=1)
        return dot/((vb*vb).sum(axis=1))