from distutils.core import setup
import sys, os.path, string, shutil

if not hasattr(sys, 'version_info') or sys.version_info < (2,2,0,'final',0):
    raise SystemExit, "Python 2.2 or later required to work with asciidata."

ver = sys.version_info
python_exec = 'python' + str(ver[0]) + '.' + str(ver[1])

def dotest():

    import Lib
    from Lib import asciifunction
    import Lib.asciidata_test
    import Lib.asciidata_numtest
    import Lib.asciidata_SExtest
    import unittest


    suite = unittest.makeSuite(Lib.asciidata_SExtest.Test_SExtractCat)
    unittest.TextTestRunner(verbosity=2).run(suite)

    suite = unittest.makeSuite(Lib.asciidata_SExtest.Test_SExtractCatII)
    unittest.TextTestRunner(verbosity=2).run(suite)

    suite = unittest.makeSuite(Lib.asciidata_SExtest.Test_SExtractCatIII)
    unittest.TextTestRunner(verbosity=2).run(suite)

    suite = unittest.makeSuite(Lib.asciidata_SExtest.Test_SExtractCatIV)
    unittest.TextTestRunner(verbosity=2).run(suite)

    suite = unittest.makeSuite(Lib.asciidata_test.Test_AsciiData)
    unittest.TextTestRunner(verbosity=2).run(suite)

    suite = unittest.makeSuite(Lib.asciidata_test.Test_AsciiDataII)
    unittest.TextTestRunner(verbosity=2).run(suite)

    suite = unittest.makeSuite(Lib.asciidata_test.Test_NullData)
    unittest.TextTestRunner(verbosity=2).run(suite)

    suite = unittest.makeSuite(Lib.asciidata_test.Test_StrangeInput)
    unittest.TextTestRunner(verbosity=2).run(suite)

    suite = unittest.makeSuite(Lib.asciidata_test.Test_SExData)
    unittest.TextTestRunner(verbosity=2).run(suite)

    suite = unittest.makeSuite(Lib.asciidata_test.Test_AsciiSort)
    unittest.TextTestRunner(verbosity=2).run(suite)

    try:
        import pyfits
        suite = unittest.makeSuite(Lib.asciidata_test.Test_AsciiFits)
        unittest.TextTestRunner(verbosity=2).run(suite)
    except:
        print 'Module pyfits is missing, skipping tests for fits.'

    try:
        import numarray
        suite = unittest.makeSuite(Lib.asciidata_numtest.Test_AsciiNumarray)
        unittest.TextTestRunner(verbosity=2).run(suite)
    except ImportError:
        print 'Module numarray is missing, skipping tests for numarray.'

    try:
        import numpy
        suite = unittest.makeSuite(Lib.asciidata_numtest.Test_AsciiNumpy)
        unittest.TextTestRunner(verbosity=2).run(suite)

        suite = unittest.makeSuite(Lib.asciidata_numtest.Test_AsciiNumpyNone)
        unittest.TextTestRunner(verbosity=2).run(suite)
    except ImportError:
        print 'Module numpy is missing, skipping tests for numpy.'


def dosetup():
    #------------------
    #  The 'normal' setup
    ret = setup(name="asciidata",
		version="1.1.1",
		description="Working with ASCII data",
		long_description="""
The asciidata package facilitates an easy working with
ASCII-tables. The ASCII data is pared and loaded into
them memory. Access to the element is provided using
matrix notation. Creting new columns, formating columns
and writing out new columns is done via high level methods.
""",
		author="Martin Kuemmel, Jonas Haase",
		author_email="mkuemmel@eso.org",
		license='LGPL',
		url="http://www.stecf.org/software/aXe/index.html",
		platforms = ["Linux","Solaris","Mac OS X"],
		packages=['asciidata',],
		package_dir={'asciidata':'Lib'})

    return ret

def main():

    if 'test' not in sys.argv:
        dosetup()
    else:
        dotest()


if __name__ == "__main__":
    main()

