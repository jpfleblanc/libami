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

	+ Sphinx >= 3.2.1 (You can install it with `pip`. Since it is python3 you may need to install `pip3`, instead of `pip`)

	+ Breathe >= 4.20.0 (You can install it with `pip` as well.)
	
3. Building:

	An example cmake command is in the `compile.sh` script.  Open the file in any editor, and change relevant file paths.
	Use a standard `CMake` procedure:

	
		::

		 $ mkdir build && cd build
		 $ cmake -DCMAKE_INSTALL_PREFIX=/path/to/install ..
		 $ make
		 $ make install


4. Creating documentation:

	Libami is documented using doxygen and Sphinx.  Documentation is set off by default and must be enabled

		::

		 $ mkdir build && cd build
		 $ cmake -DCMAKE_INSTALL_PREFIX=/path/to/install -DBUILD_DOC=ON ..
		 $ make
		 $ make install

	Documentation will then be built and installed in /path/to/install/share/doc

	Note that cmake>=3.18 is required.  As well as sphinx, the sphinx_rtd_theme and breathe python modules.  These can be ontained via

	::
		
		$pip3 install sphinx
		$pip3 install sphinx_rtd_theme
		$pip3 install breathe
		
5. Testing:

	Enabled for 'Release' builds or manually with the TEST=ON flag, are a minimal set of tests.  Simply use
		
		::

		 $ mkdir build && cd build
		 $ cmake -DCMAKE_INSTALL_PREFIX=/path/to/install -DBUILD_DOC=ON -DTEST=ON ..
		 $ make
		 $ make test
		 or
		 $ ctest

	Each test contains multiple internal checks.  If any tests should fail, more information can be gleaned by directly running the test executables in the /build/test directory.


Further information and updates will be posted on the `Github Wiki`_. 

	
------------------------
Using with your projects
------------------------

Including the header and linking to the library is required.  This can be accomplished using a standard `CMake`-based approach:

		::

		 
		  $ export libami_DIR=/install/dir/of/libami
		  $ cd /your/project/
		  $ cat >CMakeLists.txt
		  cmake_minimum_required(VERSION 3.1)
		  project(MyProject C CXX)
		  find_package(libami REQUIRED)
		  add_executable(my_program main.cpp)
		  target_link_libraries(my_program ${libami})
		  ...cont
		  $ cmake .
		  $ make




	
.. _`Github wiki`: https://github.com/jpfleblanc/libami
