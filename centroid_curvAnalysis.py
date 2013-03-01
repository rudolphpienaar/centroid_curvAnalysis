#!/usr/bin/env python

'''

    This script analyzes cumulative centroid files using the
    'C_centroid_cloud' class.
    
'''


import  os
import  sys
import  string
import  argparse
from    _common import systemMisc       as misc
from    _common import crun
from    C_centroidCloud import *

import  numpy   as np
import  pylab

import  error
import  message
import  stage
import  csv

import  fnndsc  as base
import  socket

scriptName      = os.path.basename(sys.argv[0])

class FNNDSC_CentroidCloud(base.FNNDSC):
    '''
    This class is a specialization of the FNNDSC base and geared to dyslexia
    curvature analysis.
    
    '''

    # 
    # Class member variables -- if declared here are shared
    # across all instances of this class
    #
    _dictErr = {
        'subjectSpecFail'   : {
            'action'        : 'examining command line arguments, ',
            'error'         : 'it seems that no subjects were specified.',
            'exitCode'      : 10},
        'noFreeSurferEnv'   : {
            'action'        : 'examining environment, ',
            'error'         : 'it seems that the FreeSurfer environment has not been sourced.',
            'exitCode'      : 11},
        'noStagePostConditions' : {
            'action'        : 'querying a stage for its exitCode, ',
            'error'         : 'it seems that the stage has not been specified.',
            'exitCode'      : 12},
        'subjectDirnotExist': {
            'action'        : 'examining the <subjectDirectories>, ',
            'error'         : 'the directory does not exist.',
            'exitCode'      : 13},
        'Load'              : {
            'action'        : 'attempting to pickle load object, ',
            'error'         : 'a PickleError occured.',
            'exitCode'      : 14},
    }


    def l_hemisphere(self):
        return self._l_hemi

    def l_surface(self):
        return self._l_surface

    def l_curv(self):
        return self._l_curv

    def d_centroids(self):
        return self._d_centroids

    def subj(self):
        return self._str_subj

    def surface(self):
        return self._str_surface

    def hemi(self):
        return self._str_hemi

    def curvList(self):
        return self._curvList

    def curv(self):
        return self._str_curv

    def subjDir(self):
        return "%s/%s" % (self._str_workingDir, self._str_subj)

    def analysisDir(self):
        return "%s/%s/%s/%s/%s" % \
            (self.subjDir(), _str_HBWMdir, self._str_hemi, self._str_surface, self._str_curv)

    def startDir(self):
        return self._str_workingDir

                    
    def __init__(self, **kwargs):
        '''
        Basic constructor. Checks on named input args, checks that files
        exist and creates directories.

        '''
        base.FNNDSC.__init__(self, **kwargs)

        self._lw                        = 60
        self._rw                        = 20

        # Command line arg holders
        self._str_subjectDir            = ''
        self._stageslist                = '0'
        self._hemiList                  = 'lh,rh'
        self._surfaceList               = 'smoothwm,pial'
        self._curvList                  = 'H,K'
        self._str_dataDir               = '-x'
        self._centroidTypeList          = 'pos,neg,natural'
        self._colorSpecList             = 'red,yellow,green,blue,cyan,magenta'
        self._markerSpecList            = '+,d,o,*,x,s,^'

        # List variables
        self._l_subject                 = []
        self._l_hemi                    = self._hemiList.split(',')
        self._l_surface                 = self._surfaceList.split(',')
        self._l_curv                    = self._curvList.split(',')
        self._l_type                    = self._centroidTypeList.split(',')
        self._l_color                   = self._colorSpecList.split(',')
        self._l_marker                  = self._markerSpecList.split(',')

        # Internal tracking vars
        self._str_subj                  = ''
        self._str_hemi                  = ''
        self._str_surface               = ''
        self._str_curv                  = ''
        self._str_centroidType          = ''
        self._str_markerSpec            = ''

        # Lists for tracking groups
        self._l_gidTotal                = []
        self._l_gid                     = []

        # Dictionaries for tracking data trees
        self._d_centroids               = {} # All the centroids per subject
        self._d_cloud                   = {} # This is each group's cloud
        self._d_boundary                = {} # This is each group's boundary
        self._d_poly                    = {} # This is each group's polygon boundary
        
        # Dictionaries containing all the cloud classes
        self._c_cloud                   = {}
        self._zOrderDeviation           = 3;

        self._str_workingDir            = os.getcwd()
        self._csv                       = None

        for key, value in kwargs.iteritems():
            if key == 'stages':           self._stageslist        = value
            if key == 'dataDir':
                os.chdir(value)
                self._str_dataDir       = os.path.basename(value)
            if key == 'colorSpecList':    self._l_color           = value.split(',')
            if key == 'markerSpecList':   self._l_marker          = value.split(',')
            if key == 'centroidTypeList': self._l_type            = value.split(',')
            if key == 'subjectList':      self._l_subject         = value.split(',')
            if key == 'hemiList':         self._l_hemi            = value.split(',')
            if key == 'surfaceList':      self._l_surface         = value.split(',')
            if key == 'curvList':
                self._l_curv            = value.split(',')
                self._curvList          = value

        # Read initial centroids file to determine number of subjects
        self._str_centroidFile = 'cumulative-centroids-%s.%s.%s.smoothwm.txt' % \
                (self._l_hemi[0], self._l_curv[0], self._str_dataDir)
        self._csv = csv.DictReader(file(self._str_centroidFile, "rb"), delimiter=" ", skipinitialspace=True)
        for entry in self._csv:
            self._l_subject.append(entry['Subj'])
        print self._l_subject
        
        # Build core data dictionary that contains all the centroids
        self._d_centroids               = misc.dict_init(self._l_subject)
        for subj in self._l_subject:
            self._d_centroids[subj]     = misc.dict_init(self._l_hemi)
            for hemi in self._l_hemi:
                self._d_centroids[subj][hemi]   = misc.dict_init(self._l_surface)
                for surf in self._l_surface:
                    self._d_centroids[subj][hemi][surf] = misc.dict_init(self._l_curv)
                    for curv in self._l_curv:
                        self._d_centroids[subj][hemi][surf][curv] = misc.dict_init(self._l_type)

    def centroids_read(self):
        '''
        Reads all the relevant centroid files into internal dictionary.
        '''
        for self._str_hemi in self._l_hemi:
            for self._str_surface in self._l_surface:
                for self._str_curv in self._l_curv:
                    self._str_centroidFile = 'cumulative-centroids-%s.%s.%s.%s.txt' % \
                            (self._str_hemi, self._str_curv, self._str_dataDir,
                            self._str_surface)
                    self._log('Reading centroid file: %s\n' % (self._str_centroidFile))
                    self._csv = csv.DictReader(
                                file(self._str_centroidFile, "rb"),
                                delimiter = " ",
                                skipinitialspace = True)
                    for entry in self._csv:
                        f_xn    = float(entry['xn'])
                        f_yn    = float(entry['yn'])
                        f_xp    = float(entry['xp'])
                        f_yp    = float(entry['yp'])
                        f_xc    = float(entry['xc'])
                        f_yc    = float(entry['yc'])
                        v_n     = np.array( [f_xn, f_yn] )
                        v_p     = np.array( [f_xp, f_yp] )
                        v_c     = np.array( [f_xc, f_yc] )
                        self._d_centroids[entry['Subj']][self._str_hemi][self._str_surface][self._str_curv]['neg'] = v_n
                        self._d_centroids[entry['Subj']][self._str_hemi][self._str_surface][self._str_curv]['pos'] = v_p
                        self._d_centroids[entry['Subj']][self._str_hemi][self._str_surface][self._str_curv]['natural'] = v_c
        print self._d_centroids

    def initialize(self):
        '''
        This method provides some "post-constructor" initialization. It is
        typically called after the constructor and after other class flags
        have been set (or reset).
        
        '''

        # First, this script should only be run on cluster nodes.
        lst_clusterNodes = ['rc-drno', 'rc-russia', 'rc-thunderball',
                            'rc-goldfinger', 'rc-twice']
        str_hostname    = socket.gethostname()

        # Set the stages
        self._pipeline.stages_canRun(False)
        lst_stages = list(self._stageslist)
        print self._stageslist
        for index in lst_stages:
            stage = self._pipeline.stage_get(int(index))
            stage.canRun(True)

        # Check for FS env variable
        self._log('Checking on FREESURFER_HOME', debug=9, lw=self._lw)
        if not os.environ.get('FREESURFER_HOME'):
            error.fatal(self, 'noFreeSurferEnv')
        self._log('[ ok ]\n', debug=9, rw=self._rw, syslog=False)
        
        #for str_subj in self._l_subject:
            #self._log('Checking on subjectDir <%s>' % str_subj,
                        #debug=9, lw=self._lw)
            #if os.path.isdir(str_subj):
                #self._log('[ ok ]\n', debug=9, rw=self._rw, syslog=False)
            #else:
                #self._log('[ not found ]\n', debug=9, rw=self._rw,
                            #syslog=False)
                #error.fatal(self, 'subjectDirnotExist')


    def groups_determine(self):
        '''
        Analyzes a given centroid table for all subjects and determines the
        number of groups.
        
        PRECONDITIONS
        o self._l_subject list
        
        POSTCONDITIONS
        o self._l_gidTotal
        o self._l_gid
        
        '''
        for subj in self._l_subject:
            self._l_gidTotal.append(subj[0])
        self._l_gid = sorted(set(self._l_gidTotal))
     
    def negCentroid_exists(self, str_curv):
        '''
        Returns a boolean True/False if a negative centroid exists
        for the passed str_curv.
        '''

    def clouds_gather(self):
        '''
        Populates the internal self._d_cloud dictionary across
        the input space.
        
        PRECONDITIONS
        o self._l_gid
        
        POSTCONDITIONS
        o self._d_cloud
        '''
        self._c_cloud           = misc.dict_init(self._l_gid)
        self._d_cloud           = misc.dict_init(self._l_gid)
        self._d_boundary        = misc.dict_init(self._l_gid)
        self._d_poly            = misc.dict_init(self._l_gid)
        b_firstElementPerCluster    = False
        for group in self._l_gid:
            self._d_cloud[group] = misc.dict_init(self._l_hemi)
            self._c_cloud[group] = misc.dict_init(self._l_hemi)
            self._d_boundary[group] = misc.dict_init(self._l_hemi)
            self._d_poly[group] = misc.dict_init(self._l_hemi)
            for hemi in self._l_hemi:
                self._d_cloud[group][hemi] = misc.dict_init(self._l_surface)
                self._c_cloud[group][hemi] = misc.dict_init(self._l_surface)
                self._d_boundary[group][hemi] = misc.dict_init(self._l_surface)
                self._d_poly[group][hemi] = misc.dict_init(self._l_surface)
                for surface in self._l_surface:
                    self._d_cloud[group][hemi][surface] = misc.dict_init(self._l_curv)
                    self._c_cloud[group][hemi][surface] = misc.dict_init(self._l_curv)
                    self._d_boundary[group][hemi][surface] = misc.dict_init(self._l_curv)
                    self._d_poly[group][hemi][surface] = misc.dict_init(self._l_curv)
                    for curv in self._l_curv:
                        self._d_cloud[group][hemi][surface][curv] = misc.dict_init(self._l_type)
                        self._c_cloud[group][hemi][surface][curv] = misc.dict_init(self._l_type)
                        self._d_boundary[group][hemi][surface][curv] = misc.dict_init(self._l_type)
                        self._d_poly[group][hemi][surface][curv] = misc.dict_init(self._l_type)
                        for ctype in self._l_type:
                            b_firstElementPerCluster = False
                            for subj in self._l_subject:
                                if subj[0] == group:
                                    if not b_firstElementPerCluster:
                                        self._d_cloud[group][hemi][surface][curv][ctype] = \
                                        self._d_centroids[subj][hemi][surface][curv][ctype]
                                        b_firstElementPerCluster = True
                                    else:
                                        self._d_cloud[group][hemi][surface][curv][ctype] = \
                                        np.vstack((self._d_cloud[group][hemi][surface][curv][ctype],
                                        self._d_centroids[subj][hemi][surface][curv][ctype]))
                            print self._d_cloud[group][hemi][surface][curv][ctype]
                            self._c_cloud[group][hemi][surface][curv][ctype] = \
                                C_centroidCloud(cloud=self._d_cloud[group][hemi][surface][curv][ctype])
                            self._c_cloud[group][hemi][surface][curv][ctype].confidenceBoundary_find()
                            self._d_boundary[group][hemi][surface][curv][ctype] = \
                                self._c_cloud[group][hemi][surface][curv][ctype].boundary()

    def deviation_plot(self, al_points, str_fillColor = 'red', str_edgeColor = 'black'):
        poly    = pylab.Polygon(al_points,
                            facecolor = str_fillColor,
                            edgecolor = str_edgeColor, zorder=self._zOrderDeviation)
        pylab.gca().add_patch(poly)
        return poly
        
        
    def clouds_plot(self, astr_hemi, astr_surface, astr_curv):
        '''
        Generate the actual centroid plot for given parameters
        '''
        pylab.figure()
        #pylab.axis('equal')
        pylab.grid()
        for hemi in self._l_hemi:
            for surface in self._l_surface:
                for curv in self._l_curv:
                    _d_plot     = misc.dict_init(self._l_gid)
                    for group in self._l_gid:
                        for ctype in self._l_type:
                            if ctype == 'natural': continue
                            _M_cloud = self._c_cloud[group][hemi][surface][curv][ctype].cloud()
                            _v0 = _M_cloud[:,0]
                            _v1 = _M_cloud[:,1]
                            if np.isnan(np.sum(_v0)): continue
                            _str_fileName = '%s-%s.%s.%s.%s.txt' % (group, hemi, surface, curv, ctype)
                            np.savetxt(_str_fileName, _M_cloud, fmt='%10.7f')
                            print _str_fileName
                            _d_plot[group], = plot(_v0, _v1,
                                                    color = self._l_color[int(group)-1],
                                                   marker = self._l_marker[int(group)-1],
                                                       ls = 'None',
                                                   zorder = 1)
                            self._d_poly[group][hemi][surface][curv][ctype] = \
                                self.deviation_plot(self._d_boundary[group][hemi][surface][curv][ctype],
                                               self._l_color[int(group)-1])
        pylab.show()

                
    def run(self):
        '''
        The main 'engine' of the class.

        '''
        base.FNNDSC.run(self)
            
            
def synopsis(ab_shortOnly = False):
    shortSynopsis =  '''
    SYNOPSIS

            %s                            \\
                            [--stages <stages>]                 \\
                            [-v|--verbosity <verboseLevel>]     \\
                            [--dataDir|-d <dataDir>]            \\
                            [--colorSpec|-l <colorSpec>]        \\
                            [--centroidType|-t <centroidType]   \\
                            [--hemi|-h <hemisphere>]            \\
                            [--surface|-f <surface>]            \\
                            [--curv|-c <curvType>
    ''' % scriptName
  
    description =  '''
    DESCRIPTION

        `%s' performs a centroid cloud analysis on the passed
        <Subj> <curvType> <hemi> <surface> specification.

    ARGS

        --dataDir <dataDir>
        The directory containing the centroid table files. These 
        files contain a per-subject list of centroids:

            Subj        xn      yn      xp      yp      xc      yc

        for the negative, positive, and natural centroids.

        --centroidType <centroidType>
        The "type" of centroid to analyze. One (or more) of:

                neg,pos,natural

        --colorSpec <colorSpec>
        A comma-separated string defining the colors to use for each
        centroid cloud. For example,

            'red,yellow,green,blue,cyan,magenta'

        would define six groups, with colors in order as spec'd.
    
        --hemi <hemisphere>
        The hemisphere to process. For both 'left' and 'right', pass
        'lh,rh'.

        --surface <surface>
        The surface to process. One of 'pial' or 'smoothwm'. To process both,
        use 'smoothwm,pial'.

        --curv <curvType> 
        The curvature map function stem name to analyze. The actual curvature
        file is contructed from <hemi>.<surface>.<curvType>.crv.

        --stages|-s <stages>
        The stages to execute. This is specified in a string, such as '1234'
        which would imply stages 1, 2, 3, and 4.

        The special keyword 'all' can be used to turn on all stages.


    EXAMPLES


    ''' % (scriptName)
    if ab_shortOnly:
        return shortSynopsis
    else:
        return shortSynopsis + description

def f_stageShellExitCode(**kwargs):
    '''
    A simple function that returns a conditional based on the
    exitCode of the passed stage object. It assumes global access
    to the <pipeline> object.

    **kwargs:

        obj=<stage>
        The stage to query for exitStatus.
    
    '''
    stage = None
    for key, val in kwargs.iteritems():
        if key == 'obj':                stage                   = val
    if not stage: error.fatal(pipeline, "noStagePostConditions")
    if not stage.callCount():   return True
    if stage.exitCode() == "0": return True
    else: return False


def f_blockOnScheduledJobs(**kwargs):
    '''
    A simple wrapper around a stage.blockOnShellCmd(...)
    call.
    '''
    str_blockCondition  = 'mosq listall | wc -l'
    str_blockProcess    = 'undefined'
    str_blockUntil      = "0"
    timepoll            = 10
    for key, val in kwargs.iteritems():
        if key == 'obj':                stage                   = val
        if key == 'blockCondition':     str_blockCondition      = val
        if key == 'blockUntil':         str_blockUntil          = val
        if key == 'blockProcess':
            str_blockProcess            = val
            str_blockCondition          = 'mosq listall | grep %s | wc -l' % str_blockProcess
        if key == 'timepoll':           timepoll                = val
    str_blockMsg    = '''\n
    Postconditions are still running: multiple '%s' instances
    detected in MOSIX scheduler. Blocking until all scheduled jobs are
    completed. Block interval = %s seconds.
    \n''' % (str_blockProcess, timepoll)
    str_loopMsg     = 'Waiting for scheduled jobs to complete... ' +\
                      '(hit <ctrl>-c to kill this script).    '

    stage.blockOnShellCmd(  str_blockCondition, str_blockUntil,
                            str_blockMsg, str_loopMsg, timepoll)
    return True


        
#
# entry point
#
if __name__ == "__main__":


    # always show the help if no arguments were specified
    if len( sys.argv ) == 1:
        print synopsis()
        sys.exit( 1 )

    l_subj      = []
    b_query     = False
    verbosity   = 0

    parser = argparse.ArgumentParser(description = synopsis(True))
    
    #parser.add_argument('l_subj',
                        #metavar='SUBJECT', nargs='+',
                        #help='SubjectIDs to process')
    parser.add_argument('--verbosity', '-v',
                        dest='verbosity',
                        action='store',
                        default=0,
                        help='verbosity level')
    parser.add_argument('--stages', '-s',
                        dest='stages',
                        action='store',
                        default='0',
                        help='analysis stages')
    parser.add_argument('--dataDir',
                        dest='dataDir',
                        action='store',
                        default='',
                        help='data directory containing centroid files')
    parser.add_argument('--centroidType',
                        dest='centroidType',
                        action='store',
                        default='pos',
                        help='centroid type spec to process')
    parser.add_argument('--colorSpec',
                        dest='colorSpec',
                        action='store',
                        default='',
                        help='colorSpec to process')
    parser.add_argument('--hemi', '-m',
                        dest='hemi',
                        action='store',
                        default='lh,rh',
                        help='hemisphere to process')
    parser.add_argument('--surface', '-f',
                        dest='surface',
                        action='store',
                        default='smoothwm,pial',
                        help='surface to process')
    parser.add_argument('--reset', '-r',
                        dest='b_reset',
                        action="store_true",
                        default=False)
    parser.add_argument('--curv', '-c',
                        dest='curv',
                        action='store',
                        default='H',
                        help='curvature map to process')
    args = parser.parse_args()

    OSshell = crun.crun()
    OSshell.echo(False)
    OSshell.echoStdOut(False)
    OSshell.detach(False)

    Ccloud = FNNDSC_CentroidCloud(
                        dataDir                 = args.dataDir,
                        colorSpecList           = args.colorSpec,
                        centroidTypeList        = args.centroidType,
                        stages                  = args.stages,
                        hemiList                = args.hemi,
                        surfaceList             = args.surface,
                        curvList                = args.curv,
                        logTo                   = 'CentroidCloud.log',
                        syslog                  = True,
                        logTee                  = True
                        )

    Ccloud.verbosity(args.verbosity)
    pipeline    = Ccloud.pipeline()
    pipeline.poststdout(True)
    pipeline.poststderr(True)

    stage0 = stage.Stage(
                        name            = 'CentroidCloud',
                        fatalConditions = True,
                        syslog          = True,
                        logTo           = 'CentroidCloud-process.log',
                        logTee          = True,
                        )
    def f_stage0callback(**kwargs):
        lst_subj        = []
        for key, val in kwargs.iteritems():
            if key == 'subj':   lst_subj        = val
            if key == 'obj':    stage           = val
            if key == 'pipe':   pipeline        = val
        lst_hemi        = pipeline.l_hemisphere()
        lst_surface     = pipeline.l_surface()
        lst_curv        = pipeline.l_curv()

        pipeline.centroids_read()
        pipeline.groups_determine()
        pipeline.clouds_gather()
        for pipeline._str_hemi in lst_hemi:
            for pipeline._str_surface in lst_surface:
                for pipeline._str_curv in lst_curv:
                    pipeline.clouds_plot(pipeline.hemi(), pipeline.surface(), pipeline.curv())
        os.chdir(pipeline.startDir())
        return True
    stage0.def_stage(f_stage0callback, obj=stage0, pipe=Ccloud)
    stage0.def_postconditions(f_blockOnScheduledJobs, obj=stage0,
                              blockProcess    = 'Ccloud.py')

    Ccloudlog = Ccloud.log()
    Ccloudlog('INIT: (%s) %s %s\n' % (os.getcwd(), scriptName, ' '.join(sys.argv[1:])))
    Ccloud.stage_add(stage0)
    Ccloud.initialize()

    Ccloud.run()
  
