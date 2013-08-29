#BEST OPTIMIZATION POSSIBLE I THINK. NEEDS GFORTRAN 4.7+
#env = Environment(F90='gfortran', LINK='gfortran', LINKFLAGS='', F90FLAGS='-Jobj -O3 -ffast-math -march=native -funroll-loops -fno-protect-parens -flto')

#MOST BASIC COMPILATION. SHOULD NOT BE USED BUT FOR TESTING EFFECT OF OPTIMIZATIONS
#env = Environment(F90='gfortran', LINK='gfortran', LINKFLAGS='', F90FLAGS='-Jobj')

#COMPILATION WORKING WITH OLD VERSION OF GFORTRAN
env = Environment(F90='gfortran', LINK='gfortran', LINKFLAGS='', F90FLAGS='-Jobj -O3 --fast-math -march=native -funroll-loops')


sources = Glob('src/*.f90')
objs = env.Program('meanSquaredDisplacement', sources)
