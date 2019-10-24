# Limesurf
A simulation of the mechanics of a surface as a 2D or 3D set of (visco-) elastic springs.

## Usage

```shell
    $ ./meshless config.yaml
```
In which [config.yaml](config.yaml) is an appropriate config file.

## Installation

Clone/download code in directory `<DIR>`
Then compile in `<DIR>` :

```shell
    $ cd <DIR>
    $ cmake .
    $ make 
``` 

## Config file

Config file are yaml type. See [config.yaml](config.yaml) for an example.

## Options

Pending.

## Mesh format

Format should be plain text or binary ply files. See [demo.ply](demo.ply) for an example.

## Ouput 

The program outputs ply files of the format NAME_RUN_TIME.ply, in which *NAME* is the name of the mesh provided in the config file (or the export name if given), *RUN* is the name of the run, and *T* is the time frame number.

# Serge Dmitrieff -- http://biophysics.fr
