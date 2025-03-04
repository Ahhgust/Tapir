# Virtual Tapirs!

Tapir is built for a Unix environment. That doesn't mean you can't run it on Windows, however, you just need to run Unix (on Windows).

## Virtual Box

Virtual box is a free, open source virtual machine. To unix Tapir on your Windows (or Linux) machine first:

### Install Virtual Box
This will let you use Virtual Machines on your computer..

-  Download the Virtual Box installer
   -  [External Link](https://www.virtualbox.org/wiki/Downloads)
-  Click the appropriate OS (probably Windows or Mac)
   -  Download and Install it
      -  sign the EULA and agree to everything
      -  stick to the defaults

### Configure your Windows 11 machine

Tools such as `bcl-convert` and `GLIMPSE` use special machine instructions (AVX); these instructions may not be enabled by default. AVX instructions also means that unless you're rocking an old (Intel) Mac, Tapir probably won't work.
To enable them, go to: <br><br>
`Control Panel` -> `Programs` -> `Turn Windows features on or off` (Under Programs and Features) <br><br>
And make sure that the `Virtual Machine Platform` and `Windows Hypervisor Platform` are turned off.


### Download the Tapir image
I would let this run overnight; the image is big, and because the powers that be make sharing large files exceptionally hard, I'm hosting this file on one of my home computers.
Click [here] to begin the download.

### Import the Virtual box image

-  Import the VM
   - File->Import Appliance
-  Configure the VM
   -  Adjust
      - The number of CPUs
      - The amount of RAM
      - And the disk-space
   -  To match your machine.


### Details

User:
tapir_user
Password:
4tapir

### Sanity test
When logged in as tapir_user, try:
```
~/src/Tapir/bin/bcl-convert -h
```

If you get a nice help-message, that's a good thing; if it says "Illegal instruction", that means you need to adjust how your host OS interacts with the Virtual Machine. See **Configure your Windows 11 machine** above.



