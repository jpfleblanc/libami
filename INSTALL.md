
Installation
------------

In short:
 
1. Obtaining LIBAMI:
 
	Clone Git repository:

	

	$ git clone https://github.com/jpfleblanc/libami
		
2. Prerequisites:
 
	+ C++11 compatible compiler: g++ >= 4.8 OR Intel >= 12.0 OR Clang >= 3.2

	+ CMake >= 3.18

	Optional - Documentation:

	+ Doxygen >= 1.9.0

	+ Sphinx >= 3.2.1 (You can install it by `pip`. Since it is python3 you may need to install `pip3`, instead of `pip`)

	+ Breathe >= 4.20.0 (You can install it by `pip` as well.)

3. Building:
	An example cmake command is in the `compile.sh`.  Open the file in any editor, and change relevant file paths.
	Use a standard `CMake` procedure:

		

		 $ mkdir build && cd build
		 $ sh ../compile.sh
		 $ make
		 $ make install

         
Using with your projects

Please see github page for details.

In short, use a standard `CMake`-based approach:

	

	  $ export libami_DIR=/install/dir/of/libami
	  $ cd /your/project/
	  $ mkdir build & cd build 
	  $ cmake -DCMAKE_INSTALL_PREFIX=libami_DIR ..
	  $ make



4. Creating documentation:

Libami is documented using doxygen and Sphinx.  Documentation is set off by default and must be enabled:


		$ mkdir build && cd build
		$ cmake -DCMAKE_INSTALL_PREFIX=/path/to/install -DBUILD_DOC=ON ..
		$ make
		$ make install


Documentation will then be built and installed in /path/to/install/share/doc .  Note that cmake>=3.18 is required.  As well as sphinx, the sphinx_rtd_theme and breathe python modules.  These can be ontained via:

			
		$pip3 install sphinx
		$pip3 install sphinx_rtd_theme
		$pip3 install breathe



	
	
.. _`Github wiki`: https://github.com/jpfleblanc/libami