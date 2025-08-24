import numpy as np
import pandas as pd
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
                 file_path,
                 het_nuc='N'):
        '''
        The input csv file must have columns labeled {state}w1, {state}w2 
        for states ref, A, and B

        file : str
            path to csv file with columns refw1,refw2,Aw1,Aw2,.....
        
        het_nuc : str
            Heteronucleus. Either 'N' (nitrogen) or 'C' (carbon)
        '''
        # index col needs to be the same as the ccs data and labeled 'RESI'
        try:
            df = pd.read_csv(file_path)
            if df.columns[0].upper() != 'RESI':
                print('The index column should be named "RESI".') 
            else:
                index_col = df.columns[0]
                df = pd.read_csv(file_path, index_col=index_col)
                # TODO - this isn't updating the column name  
                df.rename(columns={index_col:'RESI'},inplace=True) 
        except ValueError:
            print('Make sure that your index column is named "RESI" in all caps.')
        except FileNotFoundError:
            print("Can't find that file.")
        #print(df.head())
        self.df = df
        self.resis = df['RESI'].to_list()
        self.het_nuc = het_nuc
        if het_nuc == 'N':
            self.het_coef = 0.2
        if het_nuc == 'C':
            self.het_coef = 0.25
        # hetero nucleus shifts in column 0, H shifts in column 1 
        # TODO : add check to confirm hetero shifts are greater than H shifts
        self.states = {
        'ref': df[['refw1','refw2']].values,
        'A': df[['Aw1','Aw2']].values,
        'B': df[['Bw1','Bw2']].values
        }
        self.ref_to_B = self.get_state_distance('B')
        self.ref_to_A = self.get_state_distance('A')
        self.cos_theta = self.get_state_angle()
        self.X = self.get_fractional_activation()

    def get_state_distance(self, state):
        '''
        SI eq 5
        magnitude of vector between ref and state
        '''
        # ref to proton
        rtoH= (self.states['ref'][:,1]-self.states[state][:,1])**2
        rtohet = (self.het_coef*(self.states['ref'][:,0]-self.states[state][:,0]))**2
        return np.sqrt(rtoH+rtohet)
        
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
        shifts the protein towards activation (X > 0) or inactivation (X < 0).
        '''
        ref, A, B = list(self.states.values())
       
        # normalized vectors from ref
        vb = (B-ref)/np.linalg.norm((B-ref))
        va = (A-ref)/np.linalg.norm((A-ref))
        # row wise dot product of 2 (n_residues X 2) arrays
        dot = (va*vb).sum(axis=1)
        return dot/((vb*vb).sum(axis=1))
