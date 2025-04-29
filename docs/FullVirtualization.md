# Virtual Tapirs!

Tapir is built for a Unix environment. That doesn't mean you can't run it on Windows, however, you just need to run Unix (on Windows).

## Setting Up Virtual Box

Virtual box is a free, open source virtual machine. To unix Tapir on your Windows (or Linux) machine first:

### Minimum requirements
We recommend having
-  at least 16gb of RAM (an ABSOLUTE minimum)
-  200Gb of hard drive space
   -  and some form of solid state drive (NVME is best; traditional SSDs are fine too).
   
I was able to successfully run Tapir on my mid-grade laptop (where the minima above came from). <br>
In general, the cores you throw at some workflow, the more memory you use. If Tapir runs out of memory, use less cores.
I also made a "low memory" config (`$TAPIR/configs/config_v_2_low_mem.yaml`). This was suitable for my laptop.
e.g.,
try:
```
snakemake -c 3 (...) --configfile $TAPIR/configs/config_v_2_low_mem.yaml
```
where (...) specifies what workflow you want to run , and `-c 3` requests three CPUs.
It goes without saying that while Tapir *can* run on cheapo hardware, if you use it a fair amount, it would be wise to invest in something nicer. In general, you want more CPUs (as many as possible), more memory (~0.5Tb), and fast IO (NVMEs).
For reference, on our ~3.5 year old servers at the CHI, a miseq run can be Tapir'd in about an hour.


### Install Virtual Box
This will let you use Virtual Machines on your computer..

-  Download the Virtual Box installer
   -  [External Link](https://www.virtualbox.org/wiki/Downloads)
-  Click the appropriate OS (probably Windows or Mac)
   -  Download and Install it
      -  sign the EULA and agree to everything
      -  stick to the defaults

### Configure your Windows 11 (or Unix?) machine

Tools such as `bcl-convert` and `GLIMPSE` use special machine instructions (AVX); these instructions may not be enabled by default. AVX instructions also means that Tapir won't work on a Mac (unless you're rocking an old Intel Mac)><br>
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
		- 12Gb is a *minimum* (note this is the amount of memory dedicated to your VM. The host OS will need some too).
		- The more, the better.
      - And the disk-space
		- Aim for at least 200Gb
   -  To match your machine.

### Logging in

When you use Tapir, log in as: `tapir_user`.
```
User:
tapir_user
Password:
4tapir!
```

`tapir_user` does not have administrative privileges. `vboxuser` does. We left that password as the default: `changeme`.  Feel free to, well, change it.


### Sanity test
When logged in as tapir_user, try:
```
~/src/Tapir/bin/bcl-convert -h
```

If you get a nice help-message, that's a good thing; if it says "Illegal instruction" (and a fair amount of other text), 
that means you need to adjust how your host OS interacts with the Virtual Machine. See **Configure your Windows 11 machine** above.


### Using Tapir

When you log in as `tapir_user` the Tapir mamba environment is preloaded. 
Tapir is installed to `/home/tapir_user/src/Tapir` (i.e., that is the value of `$TAPIR`).
And remember, use the low-memory config file (tuned as appropriate).


