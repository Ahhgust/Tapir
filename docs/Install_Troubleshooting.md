# Troubleshooting

# Tapir isn't working
## 

## Glimpse
Glimpse (as provided by Tapir) requires AVX/AVX2 instructions; this may not be available on your machine (I'm looking at you, VirtualBox hypervisor). <br>
Tapir uses Glimpse v2.0.0; specifically, the release found here: <br>
https://github.com/odelaneau/glimpse/releases
<br>
Download the source-code: <br>
```
wget -O GLIMPSE.tar.gz  https://github.com/odelaneau/GLIMPSE/archive/refs/tags/v2.0.0.tar.gz
tar -xf GLIMPSE.tar.gz && rm GLIMPSE.tar.gz && cd GLIMPSE-2.0.0/
```

And follow the instructions [here](https://odelaneau.github.io/GLIMPSE/docs/installation/build_from_source)

<br>
I would add, GLIMPSE is a bit of a pain to install. Good luck.


