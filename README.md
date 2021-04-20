# DIP

The Delaunay Interpolation Program (DIP) is a C/C++ module that reads in data from the CHIPS python package and
provides interpolation on that data.

More information on CHIPS can be found here: https://github.com/Vetinar1/CHIPS

Beyond this readme file, all functions are commented and documented inside the source code itself.
Do look there if things are unclear.

## How to install DIP

There are two ways to integrate DIP into your project:
As a library (recommended if you are integrating it into a C project) or directly as part of your project
(recommended if you are integrating it into a C++ project).

### As a library

Download the source and compile it using the provided makefile. (Simply run `make`.)
This will create two shared objects libraries (*.so).
Copy the files CoolCInterface.h and CoolConst.h into your project.
Include the *.h files in your source code and use the provided functions as necessary.

At runtime, make sure that the folder containing the two *.so files is part of your `LD_LIBRARY_PATH` environment
variable so that your program can find the functions declared in the headers.


### As part of your project

Simply copy the source files into your project and use the functions directly.
For an overview over which functions these are and how to use them, see below.

No special compiler settings are required to compile DIP as part of your project.


## How to use DIP

DIP is written in C++ and makes use of C++'s object oriented programming capabilities.
It provides two classes that provide similar, but slightly different interpolation facilities.
The files CoolCInterface.* provide an adapter in case you wish to use DIP in pure C code.


### The Cool class

The Cool class, defined in CoolCool.h and CoolCool.cpp are the main class of DIP.
This class provides functions for reading DIP data, building a ball tree on that data, and running interpolation
(for more information, see thesis LINK).

Before using this class, make sure to include the CoolConst.h header file, and set the `N_MAX`, `D` and `S_MAX` constants.
`D` should be equal to the number of dimensions in your parameter space.
`N_MAX` should be equal to or slightly greater than the number of samples in your CHIPS output files.
`S_MAX` should be equal to or slightly greater than the number of simplices in your CHIPS output files.

Once these are set, instantiate a Cool object. Make sure to allocate it on the heap (using the `new`) keyword -
these objects are large and will otherwise produce a stackoverflow error.

First, set the clamping values of your parameter space using the `void set_clamp_values(double * mins, double * maxs)` method.
This function takes two array pointers as arguments, the first containing the lower limits and the second containing
the upper limits.
They automatically limit incoming interpolation values to within the defined box, to prevent errors.
Make sure to leave some padding to the edges of your sampling mesh, likewise to prevent errors (see CHIPS readme file).
**This step is important, you might otherwise get silent errors on execution.**

Afterwards, read in the CHIPS data files using the `int read_files(std::string, std::string, std::string)` method.
Its arguments are the .points, .tris and .neighbors files of the mesh that you want to read in.
This function will do a lot of calculations to set everything up for the interpolations.
For large meshes it can take several minutes to complete.

Finally, call the `void construct_btree()` method to complete the final preparation.
The Cool object is now ready for interpolations.

You can interpolate on the mesh by using the `double interpolate(double *)` method.
It takes a pointer to an array as argument. This array must contain the coordinates of the point you want to interpolate.
Make sure the coordinates are in the correct order - the same order that they appear in in the .points file.

If you wish to reuse the Cool object, you can reset it using the `void reset()` method.

If you wish to take a look at the balltree you can save it to a file using `void save_btree(std::string filename)`.
Warning: This function is not optimized and will take very long to complete for large trees.


### The CoolManager class

The CoolManager class is built on top of the Cool class and provides additional functionalities.
It is supposd to be used in the context of cosmological simulations.

CoolManager handles two Cool objects that are at different redshifts, and provides linear interpolation between the two.
It has a "high redshift" object and a "low redshift" object.
When interpolation using CoolManager, the interpolation will be passed on to these two Cool objects.
Then, the redshift dimension will be linearly interpolated between the outputs of these objects.
If the redshift moves past the interval between the two objects, the low redshift object will be pushed back and 
become the high redshift object, and a new Cool object is loaded.

To instantiate a CoolManager object a map file is needed. This file maps redshifts to sets of CHIPS data files.
Here is an example for such a map file:

```
0.0 cooldata/z0.0
0.1 cooldata/z0.1
0.2 cooldata/z0.2
0.3 cooldata/z0.3
0.4 cooldata/z0.4
0.5 cooldata/z0.5
```

It points to a set of data files called "z*" in the cooldata folder for z between 0 and 0.5.
Note the lack of file endings: Files are assumed to have the standard file endings .points, .tris and .neighbors.

Apart from that, the constructor requires values for the initial redshift interval.
The full function constructor is `CoolManager(double init_z_low, double init_z_high, std::string mapfilename)`.

After the CoolManager object is instantiated, clamping values need to be set using
`void set_clamp_values(double * mins, double * maxs)` just like for the Cool objects.
It automatically takes care of the balltrees.

Interpolation works the same as with Cool objects, too, except that redshift z needs to be given as a separate
parameter: `double interpolate(double * args, double z)`.

If z is lower than in the current interval, a new set of data is automatically loaded, as determined by the map file.
However, new "slices" can also be loaded manually using the `void push_slice(std::string)` method.
This should usually not be necessary.
If z is significantly *higher* than the current interval, an error will be thrown if CoolConst.h is configured that way.
Do not worry about rollover periods - if z is only slightly higher than z_high, CoolManager will still return
reasonable values.


## The C Interface

To use DIP with pure C code, compile it into *.so files as described above and add the CoolCInterface.h to your project.
It contains functions that enable you to use the objects described above.
Several objects can be created (how many is limited by a constant), and they can be accessed via an index that is
returned on their creation.

I recommend you take a look at CoolCInterface.cpp, it should be self-explanatory.