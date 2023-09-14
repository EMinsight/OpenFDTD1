# OpenFDTD
OpenFDTDでシミュレーションを行うためのコードです。

## 初回環境構築
### Linux開発環境
```bash
$ cd src
$ make -f Makefile_gcc clean
$ make -f Makefile_gcc　　　　OpenMP版(ofd)が作成されます
$ cd ../mpi
$ make -f Makefile_gcc clean
$ make -f Makefile_gcc　　　　MPI+OpenMP版(ofd_mpi)が作成されます
$ cd ../cuda
$ make -f Makefile_linux clean
$ make -f Makefile_linux　　　　CUDA版(ofd_cuda)が作成されます
$ cd ../cuda_mpi
$ make -f Makefile_linux clean
$ make -f Makefile_linux　　　　CUDA+MPI版(ofd_cuda_mpi)が作成されます
```

