# BLAS interfaces

Interfaces to basic linear algebra subprograms.


## Building from source

To build this project from the source code in this repository you need to have

- a Fortran compiler supporting Fortran 2008

  - GFortran 8 or newer
  - Intel Fortran 18 or newer

- One of the supported build systems

  - [meson](https://mesonbuild.com) version 0.55 or newer
  - [CMake](https://cmake.org/) version 3.9 or newer
  - [Fortran package manager (fpm)](https://github.com/fortran-lang/fpm) version 0.2.0 or newer

Get the source by cloning the repository

```
git clone https://github.com/grimme-lab/blas-interface
cd blas-interface
```

### Building with meson

To build this project with meson a build-system backend is required, *i.e.* [ninja](https://ninja-build.org) version 1.7 or newer.
Setup a build with

```
meson setup _build --prefix=/path/to/install
```

You can select the Fortran compiler by the `FC` environment variable.
To compile the project run

```
meson compile -C _build
```

Finally, you can install TOML Fortran using

```
meson install -C _build
```



### Building with CMake

While meson is the preferred way to build this project it also offers CMake support.
Configure the CMake build with

```
cmake -B _build -G Ninja -DCMAKE_INSTALL_PREFIX=/path/to/install
```

Similar to meson the compiler can be selected with the `FC` environment variable.
You can build the project using

```
cmake --build _build
```

Finally, you can install TOML Fortran using

```
cmake --install _build
```


### Building with fpm

The Fortran package manager ([fpm](https://github.com/fortran-lang/fpm)) supports the addition of this interface library as a dependency.
In the package manifest, `fpm.toml`, you can add this project dependency via:

```toml
[dependencies]
blas-interface.git = "https://github.com/grimme-lab/blas-interface"
```

Then build and test normally.

```
fpm build
fpm test
```


## License

Licensed under the Apache License, Version 2.0 (the “License”);
you may not use this file except in compliance with the License.
You may obtain a copy of the License at
http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an *“as is” basis*,
*without warranties or conditions of any kind*, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

Unless you explicitly state otherwise, any contribution intentionally
submitted for inclusion in this project by you, as defined in the
Apache-2.0 license, shall be licensed as above, without any additional
terms or conditions.
