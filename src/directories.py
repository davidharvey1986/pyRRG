import numpy as np
import os as os
def return_dirs( ):
    dirs = directories('.','.','.','.','.')
    dirs.get_dirs()
    return dirs

class directories( dict ):
    def __init__( self, data_dir, sex_files, psf_model_dir, code_dir, stilts_dir):
        self.__dict__['sex_files'] = sex_files
        self.__dict__['data_dir'] = data_dir
        self.__dict__['psf_model_dir'] = psf_model_dir
        self.__dict__['code_dir'] = code_dir
        self.__dict__['stilts_dir'] = stilts_dir

    def write_dirs( self ):
        file_obj = open('directories.cat',"wb")
        file_obj.write("DATA_DIR: %s \n" %self.data_dir)
        file_obj.write("SEX_FILES: %s \n" %self.sex_files)
        file_obj.write("PSF_MODEL_DIR: %s \n" %self.psf_model_dir)
        file_obj.write("CODE_DIR: %s \n" %self.code_dir)
        file_obj.write("STILTS_DIR: %s \n" %self.stilts_dir)

    def get_dirs( self ):
        dtypes = [('NAME', object), ('PATH', object)]
        directories = np.loadtxt('directories.cat',
                                     dtype=dtypes)
        
        self.__dict__['data_dir'] = directories['PATH'][0]
        self.__dict__['sex_files'] = directories['PATH'][1]
        self.__dict__['psf_model_dir'] = directories['PATH'][2]
        self.__dict__['code_dir'] = directories['PATH'][3]
        self.__dict__['stilts_dir'] = directories['PATH'][4]

    def check_dirs( self ):
        keys = self.keys()
        for iKey in keys:
            if not (os.path.isdir(self.__dict__[iKey])):
                raise ValueError('Cant find directory, ensure the path is correct (%s)' % self.__dict__[iKey])
            

            

    def keys(self):
        return self.__dict__.keys()

    def __getitem__(self, key): 
        return self.__dict__[key]

