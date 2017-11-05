Description
---------------
This is a repository with the standalone code of **[ichol][1]** and **[ilu][2]** functions included in GNU Octave >= 4.0.0. I developed them for the 10th edition of Google Summer of Code (2014). A set of blog posts about the development and benchmarking against Matlab's equivalents can be found [here][3].

Compiling
--------------
>**Note**: To compile the functions it is required the *mkoctfile* utility (this is included in the default octave package or in the dev package depending on your Linux distribution). 
>To build just:

    make
    
Running
------------
To launch the Octave interpreter with the functions freshly compiled use the script included:

    ./run_octave.sh

>**Note:** If you have a version of Octave >= 4.0.0 *ichol*  and *ilu* in-built versions will be override with the compiled from this repository.

 Contact
-------------
For any comment, doubt or issue write to <eduradical951@gmail.com>.

[1]: https://www.gnu.org/software/octave/doc/v4.0.0/Iterative-Techniques.html#XREFichol
[2]: https://www.gnu.org/software/octave/doc/v4.0.0/Iterative-Techniques.html#XREFilu
[3]: http://edu159-gsoc2014.blogspot.co.uk/

