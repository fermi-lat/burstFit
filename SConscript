#$Id: SConscript,v 1.3 2008/02/26 02:59:02 glastrm Exp $
Import('baseEnv')
Import('listFiles')
progEnv = baseEnv.Clone()
libEnv = baseEnv.Clone()

libEnv.Tool('burstFitLib', depsOnly = 1)
burstFitLib = libEnv.StaticLibrary('burstFit', listFiles(['src/*.cxx']))

progEnv.Tool('burstFitLib')
gtburstFitBin = progEnv.Program('gtburstFit', listFiles(['src/gtburstfit/*.cxx']))
test_burstFitBin = progEnv.Program('test_burstFit', listFiles(['src/test/*.cxx']))

progEnv.Tool('registerObjects', package = 'burstFit', libraries = [burstFitLib], includes = listFiles(['burstFit/*.h']), binaries = [gtburstFitBin],
             testApps = [test_burstFitBin], pfiles = listFiles(['pfiles/*.par']), data = listFiles(['data/*'], recursive = True))
