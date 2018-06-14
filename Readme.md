Wavemesh
========

### References ###
This code is the implementation deriving from those papers:

[VALE-04a] Valette S., Prost R., "Wavelet Based Multiresolution Analysis of Irregular Surface Meshes", IEEE Trans Visu Comp Grap, vol. 10, no. 2, pp. 113-122, 2004 .

[VALE-99a] Valette S., Kim Y S., Jung H Y., Magnin I E., Prost R., "A multiresolution wavelet scheme for irregularly subdivided 3D triangular mesh", IEEE International Conference on Image Processing ICIP'99, vol. 1, Kobe, Japan, pp. 171-174, 1999.

[VALE-04d] Valette S., Prost R., "A Wavelet-Based Progressive Compression Scheme For Triangle Meshes : Wavemesh", IEEE Trans Visu Comp Grap, vol. 10, no. 2, pp. 123-129, 2004 .

[VALE-04b] Valette S., Gouaillard A., Prost R., "Compression of 3D Triangular Meshes with Progressive Precision", Computers and Graphics, vol. 28, no. 1, pp. 35-42, 2004 .

This code is cross-platform and should compile under Linux and Window$ OS.

### License ###
This code is distributed under the GPL License
(copyright CNRS, INSA-Lyon, UCBL, INSERM.)

### Dependencies: ###
* VTK (version 6.0 or later) www.vtk.org
* CMAKE www.cmake.org

### Simple compilation HowTo under Linux ###

	git clone https://github.com/valette/Wavemesh.git
	cmake . -DCMAKE_BUILD_TYPE=Release
	make

the executables (wavemesh and others should be found under the "bin" subdirectory)

### Options ###
execute wavemesh without arguments to see the available options.

comments, suggestions : https://github.com/valette/Wavemesh/issues


