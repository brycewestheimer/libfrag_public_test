Basic Usage
============


.. code-block:: shell

    cd libfrag
    conda env create -f libfrag-dev-requirements.yml
    source activate libfrag-dev-requirements
    mkdir build
    cd build
    cmake ..
    make -j2
    make cpp-test
    make python-test
    make install
    cd examples
    ./hello_world
    cd ..
    cd benchmark
    ./benchmark_cpptools


On a windows machine this looks like:

.. code-block:: shell

    cd libfrag
    conda env create -f libfrag-dev-requirements.yml
    call activate libfrag-dev-requirements
    mkdir build
    cd build
    cmake .. -G"Visual Studio 15 2017 Win64" -DCMAKE_BUILD_TYPE=Release  ^
          -DDEPENDENCY_SEARCH_PREFIX="%CONDA_PREFIX%\Library" -DCMAKE_PREFIX_PATH="%CONDA_PREFIX%\Library"
    call activate cpptools-dev-requirements
    cmake --build . --target ALL_BUILD
    cmake --build . --target python-test
    cmake --build . --target cpp-test
    cmake --build . --target install


