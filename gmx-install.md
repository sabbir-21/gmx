# Gromacs tutorial

## Requirements
- Ubuntu 22.04.3 LTS (Download from microsoft store or dual boot) [☞Tutorial](https://youtu.be/RQKp_RA_y2k)
- [Gromacs](https://manual.gromacs.org/documentation/current/download.html) 2023.3 [☞Tutorial](https://github.com/sabbir-21/gmx/blob/main/gmx-install.md)
- [Chimera](https://www.cgl.ucsf.edu/chimera/download.html)
- [Discovery Studio](https://discover.3ds.com/discovery-studio-visualizer-download)
- [Swiss Pdb viewer](https://spdbv.unil.ch/download/binaries/SPDBV_4.10_PC.zip)
- [Gedit](```sudo apt install gedit -y```)

## Installation
### CUDA 11.8 or higher
```https://developer.nvidia.com/cuda-downloads```
CUDA 11.3, 11.5 - 11.6.1 is not compatible

For CUDA 12.3.2

```
wget https://developer.download.nvidia.com/compute/cuda/repos/ubuntu2204/x86_64/cuda-ubuntu2204.pin
```
```
sudo mv cuda-ubuntu2204.pin /etc/apt/preferences.d/cuda-repository-pin-600
```
```
wget https://developer.download.nvidia.com/compute/cuda/12.3.2/local_installers/cuda-repo-ubuntu2204-12-3-local_12.3.2-545.23.08-1_amd64.deb
```
```
sudo dpkg -i cuda-repo-ubuntu2204-12-3-local_12.3.2-545.23.08-1_amd64.deb
```
```
sudo cp /var/cuda-repo-ubuntu2204-12-3-local/cuda-*-keyring.gpg /usr/share/keyrings/
```
```
sudo apt-get update
```
```
sudo apt-get -y install cuda-toolkit-12-3
```

### GCC 10 installation
GCC 11 is not compatible
```
sudo apt update
```
```
sudo apt install gcc-10 g++-10
```
```
sudo update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-10 100
```
```
sudo update-alternatives --install /usr/bin/g++ g++ /usr/bin/g++-10 100
```
```
sudo update-alternatives --config gcc
```
```
sudo update-alternatives --config g++
```
```
gcc --version
```
```
g++ --version
```

### Gromacs Installation from source in GPU

```
wget ftp://ftp.gromacs.org/gromacs/gromacs-2023.3.tar.gz
```
```
tar xfz gromacs-2023.3.tar.gz
```
```
cd gromacs-2023.3
```
```
mkdir build
```
```
cd build
```
```
cmake .. -DGMX_BUILD_OWN_FFTW=ON -DREGRESSIONTEST_DOWNLOAD=ON -DGMX_GPU=CUDA
```
```
make
```
```
make check
```
```
sudo make install
```
```
source /usr/local/gromacs/bin/GMXRC
```

### Gromacs Installation by CONDA in CPU + gmx_MMPBSA

```
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh -b -u -p ~/miniconda3
```

Go to $Home directory /home/user/ and add ```export PATH=/home/user/miniconda3/bin:$PATH``` in the last line of .bashrc file

```
conda env create -n gmxMMPBSA --file env.yml
conda init
conda activate gmxMMPBSA
```
