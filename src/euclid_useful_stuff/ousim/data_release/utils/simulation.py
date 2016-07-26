import os
import sys
import fnmatch
import glob

import yaml
from xml.etree import ElementTree

import numpy

import pylab
from pylab import log10, cm

from IPython.core.debugger import Tracer

from astropy.io import fits
from astropy import wcs

from EuclidImageBinding import NispSurveyExposure as SurveyExposure

import catalog


class SubReleaseInfo(object):
    def __init__(self):
        self.path = None
        self.three_layer = None

    def __str__(self):
        return self.path

    @property
    def major(self):
        return self.path.split('/')[0]

    @property
    def minor(self):
        return self.path.split('/')[1]


class SubRelease(object):
    def __init__(self):
        pass
    def _add(self, release, sub_release, sub_sub_release, isInitialized):
        if hasattr(self, sub_sub_release):
            raise AttributeError('''sub-release already exists''')
        else:
            release_info = SubReleaseInfo()
            release_info.path = '%s/%s/%s' % (release,
                                              sub_release,
                                              sub_sub_release)
            release_info.three_layer = isInitialized
            setattr(self, sub_sub_release, release_info)


class _Releases(object):
    def __init__(self):
        pass
    def _add(self, release, sub_release, isInitialized=False):

        if hasattr(self, release) is False:
            setattr(self, release, SubRelease())

        if isInitialized:
            sub_sub_release = sub_release + '_3L'
        else:
            sub_sub_release = sub_release

        r = getattr(self, release)
        r._add(release, sub_release, sub_sub_release, isInitialized)
        setattr(self, release, r)

    def __iter__(self):
        """overwrite the iterator method. Returns an iterator for all the
        sub-releases"""
        releases = self
        for r in filter(lambda r: not r.startswith('_'), dir(releases)):
            sub_releases = getattr(self, r)
            for sub_r in filter(lambda r: not r.startswith('_'),
                                dir(sub_releases)):
                yield getattr(sub_releases, sub_r)

    def __str__(self):
        """overload the print method, print all the available releases."""
        return '\n'.join([release.path for release in self])

releases = _Releases()
releases._add('SR1b', sub_release='NIPf')
releases._add('SR1b', sub_release='NIPf', isInitialized=True)
releases._add('SR1c', sub_release='NIPd')
releases._add('SR1c', sub_release='NIPd', isInitialized=True)
releases._add('SC2', sub_release='NIP_test0')
releases._add('SC2', sub_release='NIP_test0', isInitialized=True)
releases._add('SC2', sub_release='NIP_test3')
releases._add('SC2', sub_release='NIP_test3', isInitialized=True)
releases._add('SC2', sub_release='NIP_R1')
releases._add('SC2', sub_release='NIP_R1', isInitialized=True)


class DetectorExposure(object):
    """Detector exposure container class that includes only the info and data
    needed by WP2500 and WP3600"""

    def __init__(self, pid, dither, detector, nirfilter):
        """constructor"""
        self.pid = pid
        self.dither = dither
        self.detector = detector
        self.nirfilter = nirfilter

        self._sci = None
        self.sci = None

        self.mask = None
        self.var = None

        self.crcount = None
        self.cr = None

        self.catalog = None
        self.planeDict = None
        self.wcs = None

        self.dataunit = None
        self.gain = None
        self.rdnoise = None
        self.saturate = None
        self.extname = None

    def set_keyword_attr_from_header_keywords(self):
        keys = ['DATAUNIT', 'RDNOISE', 'BIASLEV', 'SATURATE', 'EXTNAME',
                'GAIN']
        metadata = self._det.getMetadata()
        for key in keys:
            setattr(self, key.lower(), metadata.get(key))

    def write(self, fname):
        """"""
        pass

    def read(self, fname):
        """"""
        pass

    def bbox_wcs(self):
        """returns the corners of the detector bounding box in wcs coords.
        :return: a tuple of two arrays (ra, dec). ra and dec are the
        coordinates of corners of the bbox of the detector starting with the
        bottom left corner and moving counter-clockwise."""

        header = self._sci.header

        # the image corners (closed curve)
        pixcrd = numpy.array([# bottom left
                              [0, 0],
                              # bottom right
                              [0, header['NAXIS1']],
                              # top right
                              [header['NAXIS1'], header['NAXIS2']],
                              # top left
                              [header['NAXIS2'], 0],
                              # bottom left
                              [0, 0]], 'f4')

        ra, dec = self.wcs.wcs_pix2world(pixcrd, 1).T
        return ra, dec

    def plot_bbox(self, pid=0, dither=1, detector=0, nirfilter='Y', color='r'):
        """given the info to identify the detector uniquely plot the bounding
         box in wcs coords"""
        ra, dec = self.bbox_wcs()

        pylab.plot(ra, dec, color)

        pylab.text(ra[0], dec[0], '%s' % detector)

    def imshow(self, log=True, **kwargs):
        """show the detector science frame as an image"""

        v = self.sci.copy()
        v[ v == 0] = 1.0

        if log is True:
            v = log10(v)

        pylab.imshow(v, cmap=cm.gray, origin='lower', interpolation='nearest',
                     **kwargs)
        pylab.colorbar()


class Simulation(object):
    """Class that provides utilities and method to handle accessing pointings
    and dithers and plotting the covered reigons..etc.."""
    def __init__(self,
                 rootdir=None,
                 release=None,
                 catalog='TU',
                 three_layer=False):

        self.rootdir = rootdir
        '''The topdir containing all the simulated data'''

        self.release = release
        """The release object (identifies the release uniquely)"""

        self.nip = release.minor
        '''the sub release directory, NIPf, NIPd...etc'''

        self._cache_info_file = os.path.join(self.rootdir,
                                self.release.path,
                                'info-cache.npy')

        self.info = None
        """A summary of the important keywords collected from all the fits
        files"""

        self.catalog = catalog
        '''the dir holding the catalogs'''

        self.catmanager = None
        '''the catalog manager object'''

        self.pointings = None
        '''A pointings object'''

        self._detectors = numpy.arange(16) # asdasd
        ''' the indicies of the detectors'''

        self.pointing_ids = None
        '''An array of the unique ids of the pointings'''

        self.mdb = MDB(os.path.join(rootdir, 'MDB'))
        """the object of the MDB handler of the simulated data"""

        if os.path.isfile(self._cache_info_file):
            print('found cache info file:\n\t %s' % self._cache_info_file)
            self.info = numpy.load(self._cache_info_file)
        else:
            self.info = self.generate_info_cache(save=True)

        self.check()
        self.setup_tasks()
        # self.setup_catalog_manager()

    def find_all_data_fits_files(self):
        """returns the paths of the fits files in all the task directories. i.e
             TASK_DIR/data/*.fits
        """

        sourceDir = os.path.join(self.rootdir, self.release.path)
        lookFor = '*.fits'

        matches = []
        for root, dirnames, filenames in os.walk(sourceDir):
            for filename in fnmatch.filter(filenames, lookFor):
                if os.path.split(root)[-1] == 'data':
                    if 'CR-' not in filename:
                       matches.append(os.path.join(root, filename))
        return matches

    @staticmethod
    def set_field_id_in_info_from_RA_DEC_of_primary_header(info):
        """for each buch of RA,DEC pairs obtained from the primary header,
        assume points withing a distance of 0.1 deg to corrspond to the same
        field, thus assign a unique OBSID to them"""
        r = numpy.sqrt(info['RA']**2 + info['DEC']**2)
        from numpy import vstack,array
        from numpy.random import rand
        from scipy.cluster.vq import kmeans,vq

        data = vstack((info['RA'], info['DEC'])).T
        centroids,_ = kmeans(data, 3*3, iter=50, thresh=1e-9)
        idxes,_ = vq(data,centroids)
        for idx in idxes:
            x, y = data[idxes==idx, 0], data[idxes==idx, 1]
            pylab.plot(x, y, 'o')
            pylab.text(x.mean(), y.mean(), str(idx))
            info['OBSID'][idxes == idx] = idx
        pylab.show()

    def generate_info_cache(self, save=False):
        """from all the data fits files, collect the information from the fits
        primary header"""
        data_fits_files = self.find_all_data_fits_files()

        info = numpy.recarray(len(data_fits_files),
                              [('path', 'S1000'),
                               ('RA', 'f8'),
                               ('DEC', 'f8'),
                               ('OBSID', 'i4'),
                               ('DITHSEQ', 'i4'),
                               ('EXPNO', 'i4'),
                               ('FILTER', 'S1')])

        for i, path in enumerate(data_fits_files):
            print('.'),
            sys.stdout.flush()
            info[i]['path'] = path

            # .. todo:: fix the following bug:
            #    metadata = SurveyExposure.readFits(path).getMetadata()
            #    print metadata.keys() #! catastrophic failiur
            exposure = SurveyExposure.readFits(path)
            metadata = exposure.getMetadata()

            for name in info.dtype.names:
                if name == 'path':
                    continue
                info[i][name] = metadata.get(name)
        print()

        self.set_field_id_in_info_from_RA_DEC_of_primary_header(info)

        if save:
            numpy.save(self._cache_info_file, info)
            print('saved the cache info file to"\n\t%s' %
                  self._cache_info_file)
        return info

    @property
    def tasks_datadir(self):
        """returns the absolute path to the data directory of the pointings"""
        return os.path.join(self.rootdir, self.release.path)

    def check(self):
        """checks the values of the attributes"""
        if self.release not in releases:
            raise ValueError('''Release %s not supported...''' % self.release)

    def setup_tasks(self):
        """sets up the object for handling the available pointings (tasks)"""
        # self.pointings = PointingTasks(taskfile=self.yamlfile,
        #                                tasksdir=self.tasks_datadir)
        self.pointing_ids = numpy.unique(self.info['OBSID'])

    def setup_catalog_manager(self):
        """setup the object that handles looking up catalogs"""
        if self.release.major == 'SR1b':
            fname = os.path.join(self.rootdir, self.catalog, 'STARS', 'TU')
            self.catmanager = catalog.Manager(rootdir=fname)
            print 'setup the catalog manager'
        elif self.release.major == 'SR1c':
            print 'skipping setting up the catalog manager'

    def fetch_catalog(self, pid=0, dither=1, nirfilter='Y', cattype='STARCAT'):
        """get the catalog data the fits file data"""

        fitsfile = self.get_exposure_data_file_path(pid, dither, nirfilter)
        matchpath = os.path.abspath(
                        os.path.join(
                           os.path.join(
                              os.path.split(fitsfile)[-2],
                                      '../', '*TUStarCatalog*.fits')))
        catfilepath = glob.glob(matchpath)

        assert len(catfilepath) == 1, 'should match exactly one catalog file'

        return  catalog.Manager().catdata(catfilepath[0])

    def dither_files(self, pid=0, dither=1, nirfilter='Y'):
        """returns all the file names associated with a certain dither. i.e all
        the files mentioned in an xml file of a certain task"""
        return self.info['path'][((self.info['OBSID'] == pid) *
                                  (self.info['DITHSEQ'] == dither) *
                                  (self.info['FILTER'] == nirfilter))]

    def detector(self, pid=0, dither=1, detector=0, nirfilter='Y'):
        """give the pointing id and the dither index and the detector index
        returns a DetectorExposur object.

        :param int pid: the pointing id
        :param int dither: the index of the dither
        :param int detector: the detector index [0-15]
        :param nirfilter: not implemented
        :return: fits object of the detector and the primary HDU

        .. code-block:: python

             det = sim.detector(pid=0, dither=1, detector=0, nirfilter='Y')
             det = sim.detector(0, 1, 0, 'Y')
        """
        fname = self.get_detector_data_file_path(pid, dither, nirfilter)

        try:
            raise NotImplementedError('''not implemented''')
            retval = None
        except:

            exposure = SurveyExposure.readFits(fname)
            metadata = exposure.getMetadata()
            det = exposure.getDetector(detector)
            # sci = d.getScience().getArray()
            # Tracer()()
            # fobj = fits.open(fname)
            # primaryHUD, sci, dq = fobj[0].header,\
            #                       fobj[2*detector + 1],\
            #                       fobj[2*detector + 2]
            #
            # construct a DetectorExposure object
            retval = DetectorExposure(pid, dither, detector, nirfilter)

            # set the/ sci object and from that set copy some keywords
            retval._det = det
            retval.set_keyword_attr_from_header_keywords()

            # set the sci and mask and var layers
            retval.sci = det.getScience().getArray()
            retval.mask = det.getMask().getArray()
            retval.var = det.getVariance().getArray()

            # get the location of the CR hits (for DQ only)
            fobj = fits.open(self.get_detector_cr_data_file_path(pid, dither,
                                                                 nirfilter))
            sci_cr = fobj[2*detector + 1]
            retval.crcount = sci_cr.header['COSMICS']
            retval.cr = sci_cr.data > 0.0

            # construct the wcs object and fetch the catalog
            retval.wcs = wcs.WCS(header=sci.header)
            retval.catalog = self.fetch_catalog(pid, dither, nirfilter,
                                                'STARCAT')

        retval.planeDict = None

        return retval

    def random_detector(self):
        """returns the obejct of a random detector"""
        raise NotImplementedError("""""")

    def num_dithers_in_pointing(self, tid, nirfilter='Y'):
        """given the pointing (task) id returns the number of dithers for that
         pointing"""
        mask = (self.info['OBSID'] == tid)
        mask *= (self.info['FILTER'] == nirfilter)

        return numpy.where(mask)[0].size

    def __str__(self):
        retval = '\n'
        retval += '\nsimulation information'
        retval += '\n----------------------'
        retval += '\n%-15s: %s' % ('rootdir', self.rootdir)
        retval += '\n%-15s: %s' % ('release', self.release)
        retval += '\n%-15s: %s' % ('nip', self.nip)
        retval += '\n%-15s: %s' % ('3L', self.release.three_layer)
        retval += '\n%-15s: %s' % ('n pointings', len(self.pointing_ids))
        retval += '\n%-15s: %s' % ('pointing ids', ' '.join(map(str, self.pointing_ids)))

        retval += '\n'
        retval += '\n%-15s  %-6s %s' % ('', 'pointing_id', 'n_dithers')
        retval += '\n%-15s  %-6s %s' % ('', '           ', 'Y  J  H')
        for pid in self.pointing_ids:
            retval += '\n%-15s  %-10s' % ('', pid)
            for nirfilter in ['Y', 'J', 'H']:
                retval += '  %d' % self.num_dithers_in_pointing(pid, nirfilter=nirfilter)
        retval += '\n'
        return retval

    def get_exposure_data_file_path(self, pid=0, dither=1, nirfilter='Y'):
        """returns the path of the fits file containing the detector data"""
        return self.get_detector_data_file_path(pid, dither, nirfilter)

    def get_detector_data_file_path(self, pid=0, dither=1, nirfilter='Y'):
        """returns the path of the fits file containing the detector data"""
        fpath = self.info['path'][ (self.info['OBSID'] == pid)*
                                   (self.info['DITHSEQ'] == dither)*
                                   (self.info['FILTER'] == nirfilter)]
        if len(fpath) == 0:
            raise ValueError('''not fits file file found for
                                the specified detector''')
        else:
            assert len(fpath) == 1, 'there should be only one dectector file'
            return fpath[0]

    def get_detector_cr_data_file_path(self, pid=0, dither=1, nirfilter='Y'):
        """returns the path of the cosmic ray mask file"""
        fpath = self.get_detector_data_file_path(pid, dither, nirfilter)

        nirfilter = self.info['FILTER'][self.info['path'] == fpath][0]

        fpath_cr = fpath.replace('EUC-TEST-%s' % nirfilter,
                                 'EUC-TEST-%sCR' % nirfilter)
        assert os.path.isfile(fpath_cr)

        return fpath_cr

    def get_detector_bbox_wcs(self, pid=0, dither=1, detector=0, nirfilter='Y'):
        fpath = self.get_detector_data_file_path(pid, dither, nirfilter)
        return self._get_detector_bbox_wcs_from_file(fpath, detector=detector)

    def plot_dither_bbox_for_pointing(self, pid=0, dither=1, nirfilter='Y',
                                      color='r', detectors=[],
                                      outer_only=False):
        """give the dither index and the pointing id, plots the buonding box of
        the detectors and thus the bounding box of the dither

        :param pid: the pointing id
        :param dither: the dither index
        :param color: the colors to be cycled for the dithers
        :param detectors: the detectors to be plotted. By default, all the
        detectors are shown. An empty list indicates all detectors to be
        plotted, otherwise the detectors with the specified indicies are
        plotted.
        :param outer_only: If this is True only the outer bounding box of the
        dither is show and the individual detector bboxes are not shown
        """

        if detectors is not None:
            detectors = self._detectors if len(detectors) == 0 else detectors

        if outer_only is True:
            dither_bounds = []

        # loop over the detector (by id) and plot the bbox
        for detector in detectors:
            detector_file = self.get_detector_bbox_wcs(pid, dither,
                                                       detector, nirfilter)
            if outer_only is False:
                self._plot_detector_bbox(pid, dither, detector, nirfilter,
                                         color)
            else:
                dither_bounds.append(self.get_detector_bbox_wcs(pid,
                                                                dither,
                                                                detector,
                                                                nirfilter))

        if outer_only is True:
            dither_bounds = numpy.array(dither_bounds)

            ra = [dither_bounds[0 , 0, 0],
                  dither_bounds[3 , 0, 3],
                  dither_bounds[15, 0, 2],
                  dither_bounds[12, 0, 1],
                  dither_bounds[0 , 0, 0]]
            dec = [dither_bounds[0 , 1, 0],
                   dither_bounds[3 , 1, 3],
                   dither_bounds[15, 1, 2],
                   dither_bounds[12, 1, 1],
                   dither_bounds[0 , 1, 0]]
            pylab.plot(ra, dec, color)

    def plot_pointing_dithers(self, pid, outer_only=True, nirfilter='Y'):
        """given the ID of the pointing, plot all the dithers"""
        for dither in numpy.arange(self.num_dithers_in_pointing(pid)) + 1:
            self.plot_dither_bbox_for_pointing(pid,
                                               dither,
                                               nirfilter,
                                               outer_only=outer_only)

    def plot_pointings(self, show_all_dithers=False, outer_only=True,
                       show_sources=False, nirfilter='Y', **kwargs):
        """plots all the pointings for the simulated data.
        :param bool show_sources: if True, all the catalogs are read and the
         RA and DEC of the sources plotted.
        :param bool show_all_dithers: If True, all the dithers for each pointing
        are shown.
        :param bool outer_only: If True, the single detector boundaries are also
        shown
        """

        # get all the catalogs and plot them
        if show_sources is True:
            catalogs = self.fetch_all_catalogs(**kwargs)
            pylab.plot(catalogs.RA, catalogs.DEC, '.')

        for pid in self.pointing_ids:
            print pid
            if show_all_dithers:
                self.plot_pointing_dithers(pid,
                                           outer_only=outer_only,
                                           nirfilter=nirfilter)
            else:
                self.plot_dither_bbox_for_pointing(pid,
                                                   dither=1,
                                                   nirfilter=nirfilter,
                                                   outer_only=outer_only)

        pylab.show()

    def fetch_all_catalogs(self, cattype='STARCAT'):
        """Fetches the catalogs for all the exposures and returns them as
        a single catalog object"""

        print 'collecting all the catalogs of pointings: '
        all_catalogs = []

        # assume that all the dithers have the same catalog
        dither = 1
        for pid in self.pointing_ids:
            cat = self.fetch_catalog(pid, dither, cattype=cattype)
            all_catalogs.append(cat)
        all_catalogs = numpy.hstack(all_catalogs)
        return numpy.recarray(all_catalogs.shape,
                              dtype=all_catalogs.dtype,
                              buf=all_catalogs.data)

class MDB(object):
    """handle accessing the files in the MDB directory

    .. code-block:: python

        mdb = MDB(rootdir='~/data/simulated/MDB')
    """
    def __init__(self, rootdir):
        """initialize the MDB object.
        :param rootdir: The path of the MDB directory.
        """
        self._rootdir = None
        """backing variable for the rootdir property"""

        self.rootdir = rootdir

    @property
    def rootdir(self):
        return self._rootdir

    @rootdir.setter
    def rootdir(self, p):
        """setter for the _rootdir backing attribute"""
        self._rootdir = os.path.expanduser(p)
        assert os.path.exists(self._rootdir)

    def distortions(self, coordinate=None, band=None):
        """return the content of the fits file of the specific coordinate and
        the specific filter
        :param coordinate: The coordinate string 'x' or 'y'
        :param band: The filter H, J, or Y
        """

        assert coordinate in ['x', 'y']
        assert band in ['H', 'J', 'Y']

        distortion_file_path = os.path.join(self.rootdir,
                                            'distortion_%s_%s.fits' % (
                                                coordinate,
                                                band))

        assert os.path.exists(distortion_file_path)

        f = fits.open(distortion_file_path)

        print 'opened distortion file:\n\t%s' % distortion_file_path

        return f[0]
