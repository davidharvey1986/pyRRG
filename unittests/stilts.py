import match_cat
import pyfits as fits

def test_all( ):
    '''
    Test the stilts path and make sure i am able to
    match catalogues
    '''

    path_test()
    catlalogue_test()


def path_test():
    '''
    Test sextractor is in the path
    '''
    try:
        sex_dir = '/'.join(subprocess.check_output(['which','stilts.sh']).split('/')[:-1] )
    except:
        raise ImportError("Path Test Failed: Cant find SExtractor in path")
    
