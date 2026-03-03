import cluster
import plots
import svd
import pandas as pd

# This is intended to give the option to run chesca interactively.
# output the dendrogram, svd, chespa (if reference state is available), and input parameters used
# provide option to list the cluster ids you want to remove and run again (with new parameters potentially)
# if pymol (try: pymol -r) and structure is available, automatically generate the clusters
# in the viewer.

run = True

while run == True:
    file_path = input('Enter the path to your combined chemical shift file:  ')
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
    print(df.head())
    # SVD object
    dims = svd.SVD(df)
    cutoff = input("Enter the correlation coefficient cutoff for agglomerative clustering:  ")
    