An educative tool to visually inspect concepts such as temperature, thermal equilibrium and the probability density distribution of velocities in a gas. Released under GNU GPL 3 by Jose Luis Blanco (University of Málaga).

The basics of the theory behind the simulator was explained, along with a few videos, in these (Spanish) popular Science blogs:
  * http://www.ciencia-explicada.com/2012/06/que-es-la-temperatura-una-explicacion.html
  * http://amazings.es/2012/06/26/que-es-la-temperatura-una-explicacion-traves-de-animaciones/


## 1. Usage ##

These keys are recognized as commands during simulations:

```
 0   : Add new static particle.
 1-4 : Launch particles of increasing speeds.
 6   : Launch a very fast (x100) particle.
 8-9 : Launch a lot of new particles.
 t   : Toogle tracking of random-walk.
 x/z : Simulation speed up/down (limited by CPU time).
 +/- : Increase/decrease point particle size.
```

Sample screenshot:

![http://maxwell-boltzmann-simulator.googlecode.com/svn/trunk/screenshot1.png](http://maxwell-boltzmann-simulator.googlecode.com/svn/trunk/screenshot1.png)


## 2. Compiling: Ubuntu ##

  * Install CMake and [MRPT](http://www.mrpt.org/MRPT_in_GNU/Linux_repositories).
  * cmake . && make
  * Launch 2D or 3D simulators with: "./boltz-2d" or "./boltz-3d"


## 3. Compiling: Windows ##

  * Install CMake and MRPT.
  * Generate Visual Studio projects from cmake-gui, compile from the IDE.