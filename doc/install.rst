============
Installation
============

For detailed instructions, please see `Github Wiki`_

In short:
 
1. Obtaining LIBAMI:
 
	Clone Git repository:

	::

	$ git clone https://github.com/jpfleblanc/leblanc_codes
		
2. Prerequisites:
 
	+ C++11 compatible compiler: g++ >= 4.8 OR Intel >= 12.0 OR Clang >= 3.2

	+ CMake >= 3.18

	Optional:

	+ Doxygen >= 1.9.0

	+ Sphinx >= 3.2.1 (You can install it by `pip`. Since it is python3 you may need to install `pip3`, instead of `pip`)

	+ Breathe >= 4.20.0 (You can install it by `pip` as well.)

3. Building:
	An example cmake command is in the `compile.sh`.  Open the file in any editor, and change relevant file paths.
	Use a standard `CMake` procedure:

		::

		 $ mkdir build && cd build
		 $ sh ../compile.sh
		 $ make
		 $ make install

         
Using with your projects

Please see github page for details.

In short, use a standard `CMake`-based approach:

	::

	  $ export libami_DIR=/install/dir/of/libami
	  $ cd /your/project/
	  $ cat >CMakeLists.txt
	  cmake_minimum_required(VERSION 3.1)
	  project(MyProject C CXX)
	  find_package(libami REQUIRED)
	  add_executable(my_program main.cpp)
	  target_link_libraries(my_program ${libami})
	  ^D
	  $ cmake .
	  $ make



4. Creating documentation:

Libami is documented using doxygen and Sphinx.  Documentation can be built using 

	::
	
		$make LibAMI-docs
		$make Sphinx 



Further information and updates will be posted on the `Github Wiki`_. 

	
	
.. _`Github wiki`: https://github.com/jpfleblanc/leblanc_codes
