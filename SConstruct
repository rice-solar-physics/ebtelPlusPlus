"""
Build file for ebtel++
"""
import sys
import os

AddOption('--debug_compile', dest='debug_compile', action='store_true', help='Turn off optimizations for debugging.')
# Allow user to specify custom path to libs and headers
# If all else fails, just pass the absolute paths with the libs flag
AddOption('--libpath', dest='libpath',type='string',nargs=1,action='store',default=None,
          help='Comma-separated list of custom library paths if the defaults do not work.')
AddOption('--libs', dest='libs', type='string', nargs=1, action='store', default=None,
          help='Comma-separated list of libraries. This is most useful for modifying the name of the Boost library.')
AddOption('--includepath', dest='includepath', type='string', nargs=1, action='store', default=None,
          help='Comma-separated list of custom include paths if defaults do not work.')
AddOption('--linkflags', dest='linkflags', type='string', nargs=1, action='store', default=None,
          help='General options passed to the linker. Use this to pass -rpath to the linker if needed.')


subdirs = ['source']
cxx_flags = ['-std=c++14',]
try:
    CXX = os.environ['CXX']
except KeyError:
    CXX = 'g++'
if GetOption('debug_compile'):
    cxx_flags += ['-g', '-Wall',]
else:
    cxx_flags += ['-O3']
env = Environment(CXX=CXX, CXXFLAGS=cxx_flags)

if 'darwin' in sys.platform:
    print("Using Mac OS X compile options.")
    env.Append(CPPPATH=['/opt/local/include', '/usr/include/malloc'])
    env.Append(LIBS=['boost_program_options-mt'])
    env.Append(LIBPATH=['/opt/local/lib'])
elif 'linux' in sys.platform:
    print("Using Linux compile options.")
    env.Append(CPPPATH=['/usr/include'])
    env.Append(LIBS=['m', 'stdc++', 'boost_program_options'])
    env.Append(LIBPATH=['/usr/lib/x86_64-linux-gnu'])
else:
    print("Unrecognized platform. Using Windows compile options.")
    env.Append(CPPPATH=['/usr/local/include'])
    env.Append(LIBS=['boost_program_options'])
    env.Append(LIBPATH=['/usr/local/lib'])

if GetOption('libpath'):
    env['LIBPATH'] = GetOption('libpath').split(',')
if GetOption('includepath'):
    env['CPPPATH'] = GetOption('includepath').split(',')
if GetOption('libs'):
    env['LIBS'] = GetOption('libs').split(',')
if GetOption('linkflags'):
    env['LINKFLAGS'] = GetOption('linkflags')

allobjs = []
for sd in subdirs:
    consfile = os.path.join(sd, 'SConscript')
    allobjs += env.SConscript(consfile, exports=['env'])

# TODO will this path always be resolved correctly?
if not os.path.exists('bin'):
    os.makedirs('bin')

env.Program('bin/ebtel++.run', allobjs)
