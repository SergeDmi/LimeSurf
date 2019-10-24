# Limesurf
A simulation of the mechanics of a surface as a 2D or 3D set of (visco-) elastic springs.

## Usage

```shell
    $ ./meshless config.yaml
```
In which [config.yaml](config.yaml) is an appropriate config file.

## Installation

Clone/download code in directory `<DIR>` :

 ```shell 
    $ git checkout https://github.com/SergeDmi/LimeSurf.git <DIR>
    $ cd <DIR>
    $ git submodule update --init --recursive
```

Then compile in `<DIR>` :
```shell
    $ cd <DIR>
    $ cmake .
    $ make 
``` 

For this you might need to install make, cmake, and the boost libraries. To do so on Ubuntu, use :
```shell
    $ sudo apt-get install make cmake
    $ sudo apt-get install libboost-all-dev
```
On a Mac, you can use *brew install* if you are using [Homebrew](https://brew.sh) as your package manager.
On Windows, you can install Ubuntu on a virtual machine, see this [tutorial](https://brb.nci.nih.gov/seqtools/installUbuntu.html) for instance.

## Config file

Config file are yaml type. See [config.yaml](config.yaml) for an example.

## Options

Pending.

## Mesh format

Format should be plain text or binary ply files. See [demo.ply](demo.ply) for an example.

## Ouput 

The program outputs ply files of the format NAME_RUN_TIME.ply, in which *NAME* is the name of the mesh provided in the config file (or the export name if given), *RUN* is the name of the run, and *T* is the time frame number.

# Serge Dmitrieff -- http://biophysics.fr
