import os
import fnmatch
from astropy.io import fits
from IPython.core.debugger import Tracer

class Manager(object):
    def __init__(self, rootdir=None):
        """
        .. code-block:: python

             m = Manager(rootdir='../SIMCHallenge1b/stars/TU')
             fname = m.find(4938, 'STARCAT')
             fname = m.find(4938, 'STARSPC')
        """
        self.rootdir = rootdir

    def taskid_type_from_fname(self, fname):
        """returns the taskid and the catalog type from the filename"""
        basenameSplit = os.path.split(fname)[-1].split('-')
        cattype, taskid = basenameSplit[2], basenameSplit[3]
        return cattype, int(taskid)

    def find(self, taskid, cattype):
        """return the full path of the file that matches the taskid and the
        cattype. """
        matches = []

        for root, dirnames, filenames in os.walk(self.rootdir):
          for filename in fnmatch.filter(filenames, '*.fits'):
              typeGuess, tidGuess = self.taskid_type_from_fname(filename)
              if cattype in typeGuess and taskid == tidGuess:
                  matches.append(os.path.join(root, filename))

        if len(matches) > 1:
            raise ValueError('''It seems the catalog file found is not
            unqiue. More than one catalog found:\n\t\t%s
                            ''' % '\n\t\t'.join(matches))
        if len(matches) == 0:
            raise ValueError('''the specified catalog could not be located
            in the directory:\n\t%_s''' % self.rootdir)

        return matches

    @staticmethod
    def catdata(fname):
        """return the content of the catalog given the path"""
        return Catalog(fname).read().data

    def fetch(self, taskid, cattype):
        """returns the catalog data of the specified task given the catalog
        type"""
        fname = self.find(taskid, cattype)[0]
        print 'catalog filename:\n\t%s' % fname
        return self.catdata(fname)

class Catalog(object):

    def __init__(self, fname):
        self.fname = fname

    def read(self):
        hdulist = fits.open(self.fname)
        primaryHDU, binTableHdu, = hdulist
        return binTableHdu
