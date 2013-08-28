# meanSquaredDisplacement

## Purpose

Computes the mean squared displacement of particles along a trajectory.
It takes the periodic boundary counditions into account.

## Background

The purpose of this program is to compute the mean squared displacement (MSD) of a series of equivalent particles that have positions
evolving in time.

The MSD we calculate is defined as  
<a href="http://imgur.com/2mvwmPO"><img src="http://i.imgur.com/2mvwmPO.jpg?1" title="Hosted by imgur.com" /></a>

where **q** is the three-dimension position vector. It is averaged over all sites, *i*, and origins of time, *t0*.

In a nutshell, the MSD of a time *t* is the ensemble average of the squared Euclidean distance between two positions separated by time *t*

## Author

Written by Maximilien Levesque, while in postdoc in the group of Mathieu Salanne at  
UPMC Univ Paris 06, CNRS, ESPCI, UMR 7195, PECSA, F-75005 Paris, France

## Thanks

* Marie Jardat, UPMC, PECSA, Paris, France, for discussions, extensive testing and bug reports of the beta (unshared) versions.  
* Xudong Zhao, UPMC, PECSA, Paris, France, for providing test cases and for discussions.

## How to compile it

You need to install `scons`, which is a newer and better and easier equivalent to the much complicated `make`.
Please visit [the scons website](www.scons.org) for download, or you should better use the repositories of your own distribution:  
Under Ubuntu: `sudo apt-get install scons`  
Under Fedora: `sudo yum install scons`  
and so on for other distros.

Once you're in the directory where you downloaded `meanSquaredDisplacement`, just type  
```
$ scons
```  
In a nutshell, scons will look at the file called `SConstruct` I built for you and compile everything smartly.

Of course, you may compile the simple fortran file(s) by yourself if you prefer the complicated ways.

## How to use meanSquaredDisplacement

### Inputs
All you need is a file, whatever its name, containing all your positions in a format  
``` 
x(1,t) y(1,t) z(1,t)  
x(2,t) y(2,t) z(2,t)  
x(i,t) y(i,t) z(i,t)  
...    
x(Nat,t) y(Nat,t) z(Nat,t)    
x(1,t+1) y(1,t+1) z(1,t+1)  
x(2,t+1) y(2,t+1) z(2,t+1)  
x(i,t+1) y(i,t+1) z(i,t+1)  
...  
x(Nat,t+1) y(Nat,t+1) z(Nat,t+1)  
...  
...  
x(Nat,t+Nstep) y(Nat,t+Nstep) z(Nat,t+Nstep)  
```  

where `Nat` is the total number of atoms you have in your supercell,  
and `Nstep` is the number of timesteps in your trajectory.

In other words, you print the coordinates of all sites for a given timestep, then for the next one, etc.  
You do not print blank lines anywhere. You can add as many spaces you want between the columns.

At the end, you *must* have 3 columns and `Nat x Nstep` rows.

That's all you need.

Note that if you have one-dimensional (two-dimensional) positions at each time step, you'll have to add 2 (1) column filled with zero.

## Execution

The executable is waiting for arguments:  
1. `Lz`, the length of the supercell in *x* direction, in the same units as the positions  
2. `Ly`  
3. `Lz`  
4. `Nat`, defined above  
5. `filename` of the trajectory in the format discussed above.  
  
So, you have to execute:  
`$ meanSquaredDisplacement 56.0 37.43 288.1 800 ./analysis/MSD/positions.out`

## Outputs

Execution of `meanSquaredDisplacement` will result in a single file: `msd.out`.  
This file has a very simple ASCII format:
``` 
1 1.12931233  
2 1.24034134  
3 1.5908O123  
...  
Nstep 123987.12383  
``` 
where the first column is the timestep, and the second column is the mean squared displacement.  

### Units
Both timesteps and mean squared displacement are in units of your simulation.  
If your trajectory file contains milliseconds (ms) and (Angstroms), `msd.out` will be in `ms Angstroms²`.


## Disclaimer

This program has been thoroughly tested and validated.
I would not share a tool that is not production ready. I am very confident in this program.
*Nevertheless*, this program may contain bugs or restrictions that would lead to unexpected results.
Please, read carefully this readme file, and test it thoroughly on data for which you already know the results.

Comments, bug reports or even thanks ;) are welcome!
