# -*- python -*-
# $Id: SConscript,v 1.8 2009/07/20 20:13:14 peachey Exp $
# Authors: James Peachey <James.Peachey-1@nasa.gov>
# Version: burstFit-02-02-03
Import('baseEnv')
Import('listFiles')
progEnv = baseEnv.Clone()
libEnv = baseEnv.Clone()

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

