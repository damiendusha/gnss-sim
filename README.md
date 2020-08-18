# GNSS-SIM

This is a branch of the classic GPS simulation program 
[gps-sdr-sim](https://github.com/osqzss/gps-sdr-sim). The key differences
between this and the original are:

 - Split out the single source into multiple C++ files; and
 - Added a facility for constant per-satellite signal strength.
 
The latter capability is to facilitate testing with known GPS waveforms. There
is currently a hack so that this overrides the "proper" link budget calculation.

Of course, many features are planned and very few will make it to fruition!

### Building

```
cmake .
make
```

### Tested Platforms

Tested on Ubuntu 20.04 generated with an Ettus N200 SDR to a U-Blox receiver.
 
 
### License

Copyright &copy; 2020 Damien Dusha__
Distributed under the [MIT License](http://www.opensource.org/licenses/mit-license.php).

Derived from (gps-sdr-sim](https://github.com/osqzss/gps-sdr-sim):

Copyright &copy; 2015-2020 Takuji Ebinuma
Distributed under the [MIT License](http://www.opensource.org/licenses/mit-license.php).
