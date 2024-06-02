# rCTCBN
An adaptation of the CT-CBN program for R
## Installation
`rCTCBN` can be installed from GitHub using
```R
# install.packages("devtools")
devtools::install_github("rockwillck/rCTCBN")
```
Run `?ctcbn` to learn how to use the function.
## Dependencies
### For Mac
Because Mac's distribution of `clang` does not support OpenMP out of the box, the following shell command is necessary.  
```console
foo@bar:~$ curl -O https://mac.r-project.org/openmp/openmp-16.0.4-darwin20-Release.tar.gz
foo@bar:~$ sudo tar fvxz openmp-16.0.4-darwin20-Release.tar.gz -C /
```
## Testing
- [x] M1 Mac
- [x] M2 Mac
- [ ] Linux
- [ ] Windows
## Roadmap
- [x] Base method
- [ ] Parameters
- [ ] Documentation
