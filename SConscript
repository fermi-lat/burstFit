# -*- python -*-
# $Id: SConscript,v 1.12 2015/03/21 05:38:57 jchiang Exp $
# Authors: James Peachey <James.Peachey-1@nasa.gov>, Joe Asercion <joseph.a.asercion@nasa.gov>
# Version: burstFit-02-03-01
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

