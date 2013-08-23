env = Environment(F90='gfortran', LINK='gfortran', LINKFLAGS='', F90FLAGS='-Jobj -O3 -ffast-math -march=native -funroll-loops -fno-protect-parens -flto')
#env = Environment(F90='gfortran', LINK='gfortran', LINKFLAGS='', F90FLAGS='-Jobj -Og')
sources = Glob('src/*.f90')

# The next line is the actual code that links the executable. env.Program generates an executable
objs = env.Program('meanSquaredDisplacement', sources)
