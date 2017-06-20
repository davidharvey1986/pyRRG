import subprocess

def check_external_packages():

    try:
        stilts_path = subprocess.check_output(['which','stilts.sh'])
    except:
        raise ValueError('Cannot find STILTS please install and ensure it is in the shell path')

    
    try:
        stilts_path = subprocess.check_output(['which','sex'])
    except:
        raise ImportError('Cannot find SExtractir please install and ensure it can be called with "sex"')

    try:
        import pickle as pkl
    except:
        raise ImportError('Cannot find pickle, plesae install')
