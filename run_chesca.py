import cluster
import plots
import svd
import pandas as pd

# This is intended to give the option to run chesca interactively.

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
    