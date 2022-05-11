'''
Utility Functions to set up FDASCC Package
'''

def set_r_home(path_to_R='C:\Program Files\R\R-4.2.0'):
    import os
    if 'r_home' in os.environ:
        print(f"r_home is already set to {os.environ['r_home']}")
    else:
        os.environ['r_home'] = path_to_R
        print(f"r_home environment variable set to {path_to_R}")


def install_fdascc(path_to_fdascc="C:/Users/alley/Desktop/folders/courses-and-certifications/phd-cs/research/1d_fda/FDASCC_py/FDASCC/FDASCC"):
    import rpy2
    from rpy2 import robjects
    import rpy2.robjects.packages as rpackages
    
    utils = rpackages.importr('utils')
    # utils.chooseCRANmirror(ind=1)
    utils.install_packages(path_to_fdascc, 
                           repos=rpy2.rinterface.NULL,
                           type="source")
    print("FDASCC successfully installed")