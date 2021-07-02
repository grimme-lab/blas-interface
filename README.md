# BLAS interfaces

Interfaces to basic linear algebra subprograms.


## Usage

This project supports the Fortran package manager ([fpm](https://github.com/fortran-lang/fpm)) as build system.
To use this interfaces in your project add `blas-interface` as dependency to your package manifest.

```toml
[dependencies]
blas-interface.git = "https://github.com/awvwgk/blas-interface"
```

The `blas` module and the `blas_pure` module become importable.
The latter declares all BLAS procedures with explicit intent to allow using them in a `pure` context, even if the procedures themselves are not declared as `pure`.


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
