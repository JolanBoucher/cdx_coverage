# cdx_coverage

`cdx_coverage` computes sequencing coverage over a linearized pangenome
representation using:

- a custom CDX index; and
- read alignments stored in VG's GAM format.

The GAM file must have been generated from the same pangenome graph used to
build the CDX index. In particular, the graph node identifiers referenced by
the GAM alignments must match the node identifiers stored in the CDX index.

## Requirements
#### Primary dependency
Requirements:
- Python >= 3.10
- numpy
- matplotlib

- [libvgio](https://github.com/vgteam/libvgio)

  Used to deserialize GAM files and access VG Protobuf alignment types such as
  `vg::Alignment`, `vg::Mapping`, `vg::Position`, and `vg::Edit`.

#### libvgio dependencies

- [Protobuf](https://github.com/protocolbuffers/protobuf)
- [HTSlib](https://github.com/samtools/htslib), version 1.10 or newer
- POSIX threads
- OpenMP

#### Build tools

- A C++17-compatible compiler
- [cmake](https://cmake.org/), version 3.16 or newer
- [pkg-config](https://wwworg/wiki/Software/pkg-config/)

## Installing dependencies
### macOS (Silicon ARM64)

Install the dependencies using Homebrew:

```bash
brew install cmake protobuf htslib libomp pkg-config
```

### Ubuntu 20.04 or newer

Install the compiler, CMake, Protocol Buffers, HTSlib, OpenMP, and the other
required development packages from the Ubuntu repositories:

```bash
sudo apt update
sudo apt install -y \
    build-essential \
    cmake \
    git \
    pkg-config \
    protobuf-compiler \
    libprotobuf-dev \
    libhts-dev \
    libomp-dev \
    zlib1g-dev 
```