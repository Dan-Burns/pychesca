# functions to prepare pymol files that let you visualize the chesca results
# on the structure.


def clusters_to_pymol(df,output):
    '''
    take a df of residue indices and cluster labels
    and make pymol selections
    '''
    colors = ['red','blue','yellow','black','white','purpleblue','magenta','orange','neptunium','zinc','radium','pink',
             'purple','lighblue','teal','marine','brightorange','chocolate','sand','slate']
    
    clusters = {}
    # Get a set of the cluster ids {1,2,3,...}
    for cluster_label in set(df['cluster']):
        # Add the residue ids to the dictionary that correspond to the cluster id
        clusters[int(cluster_label)] = list(df.loc[df['cluster'] == cluster_label].index)
    
    with open(output,'w') as f:
        f.write('set sphere_scale, 0.8\n')
        for i, cluster in enumerate(clusters):
            sel = ''
            for resid in clusters[cluster]:
                # int(float) to catch Aayushi's residue decimal naming scheme to differentiate CH/NH/ILE methyls etc.
                sel += f'{int(float(resid))}+'
                #sel += f'{(resid)}+'

            f.write(f'select c_{cluster}, resi {sel}\n') 
            f.write(f'color {colors[i]}, c_{cluster}\n show spheres, c_{cluster}\n')
            

# CHESCA clusters are often in agreement with local RMSD maxima -
# include RMSD/ CHESCA visualization

# also possible for no change in RMSD + CHESCA 
# or no CHESCA but large changes
# look for changes in rates
