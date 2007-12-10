import glob,os,platform

Import('baseEnv')
Import('listFiles')
progEnv = baseEnv.Clone()
libEnv = baseEnv.Clone()

burstFitLib = libEnv.StaticLibrary('burstFit', listFiles(['src/*.cxx']))

progEnv.Tool('burstFitLib')
gtburstFitBin = progEnv.Program('gtburstFit', listFiles(['src/gtburstfit/*.cxx']))

progEnv.Tool('registerObjects', package = 'burstFit', libraries = [burstFitLib], includes = listFiles(['burstFit/*.h']), binaries = [gtburstFitBin], pfiles = listFiles(['pfiles/*.par']))