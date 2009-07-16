# -*- python -*-
# $Id: SConscript,v 1.6 2009/07/15 23:11:28 jrb Exp $
# Authors: James Peachey <James.Peachey-1@nasa.gov>
# Version: burstFit-02-02-02
Import('baseEnv')
Import('listFiles')
progEnv = baseEnv.Clone()
libEnv = baseEnv.Clone()

libEnv.Tool('burstFitLib', depsOnly = 1)
burstFitLib = libEnv.StaticLibrary('burstFit', listFiles(['src/*.cxx']))

progEnv.Tool('burstFitLib')
gtburstFitBin = progEnv.Program('gtburstFit', listFiles(['src/gtburstfit/*.cxx']))
test_burstFitBin = progEnv.Program('test_burstFit', listFiles(['src/test/*.cxx']))

progEnv.Tool('registerTargets', package = 'burstFit',
             staticLibraryCxts = [[burstFitLib,libEnv]],
             includes = listFiles(['burstFit/*.h']),
             binaryCxts = [[gtburstFitBin, progEnv]],
             testAppCxts = [[test_burstFitBin,progEnv]],
             pfiles = listFiles(['pfiles/*.par']),
             data = listFiles(['data/*'], recursive = True))

