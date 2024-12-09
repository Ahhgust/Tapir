#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;

my $sep = "\t";
my $basename;
my $post = 0;
my $echo = "";

my $header=0;

GetOptions('h' => \$header,
    's=s' => \$sep, # separator (tabs)
           'e=s' => \$echo, # echos additional (constant) information
	   'p' => \$post, # prints filename before (pre, default) or post
	   'b:s' => \$basename); # prints just the basename of the file (not the path). behavior as per the unix utility basename...


my $prevARGV='';
my $fn=''; # filename

my $printFirst=$header;

while (<>){
    
    if (defined $basename) {
        if ($ARGV ne $prevARGV) {
            my @s = split /\//, $ARGV;
            $fn = $s[-1]; # get the filename w/o the path (unix basename style)
            
            if ($basename ne '') { # optionally trim off a file extension
                $fn =~ s/$basename$//;
            }
            
            $prevARGV = $ARGV;
            next if $header;
        }
    } else {
        $fn = $ARGV;
    }

    # used to print ammend the header...
    if ($printFirst) {
        $fn = "Filename";
        $printFirst=0;
        $prevARGV = $ARGV;
    } elsif ($header && $prevARGV ne $ARGV) {
        $prevARGV = $ARGV;
        next;
    }    
    
    if ($post) {
        chomp;
        if ($echo ne '') {
            print $_ , $sep , $fn ,  $sep , $echo ,"\n";
        } else {
            print $_ , $sep , $fn , "\n";
        }
    } else {
        if ($echo ne '') {
            print $echo , $sep , $fn , $sep , $_;
        } else {
            print $fn , $sep , $_;
        }
    }
}

