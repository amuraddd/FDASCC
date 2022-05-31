'''
Utility Functions to set up FDASCC Package
'''
def setup_fdascc(path_to_R='C:/Users/alley/AppData/Local/R-3.6.3', path_to_fdascc="C:/Users/alley/Desktop/folders/courses-and-certifications/phd-cs/research/FDASCC/FDASCC/FDASCC"):
    import os
    if 'r_home' in os.environ:
        print(f"r_home is already set to {os.environ['r_home']}")
    else:
        os.environ['r_home'] = path_to_R
        print(f"r_home environment variable set to {path_to_R}")
        
        
    import rpy2
    from rpy2 import robjects
    import rpy2.robjects.packages as rpackages
    from rpy2.robjects.vectors import StrVector # R vector of strings
    
    utils = rpackages.importr('utils')
    utils.chooseCRANmirror(ind=1)
    # install packages
    packnames = ('gtools', 'splines2', 'pracma', 'RSpectra')

    # install only packages which need to be installed
    names_to_install = [x for x in packnames if not rpackages.isinstalled(x)]
    if len(names_to_install) > 0:
        utils.install_packages(StrVector(names_to_install))
    
    # install FDASCC
    if rpackages.isinstalled('FDASCC'):
        print("FDASCC Installed")
    else:
        utils.install_packages(path_to_fdascc,
                               # lib=path_to_R+'\library',
                               repos=rpy2.rinterface.NULL,
                               type="source")
        print("FDASCC Successfully Installed")