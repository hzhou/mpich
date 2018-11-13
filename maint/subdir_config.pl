#!/usr/bin/perl
use strict;
our @dir_list;
our %opts=(j=>1);

foreach my $a (@ARGV){
    if(-d $a){
        push @dir_list, $a;
    }
    elsif($a=~/-j(\d+)/){
        $opts{j}=$1;
    }
}
if(-d ".config"){
    system "rm -r .config";
}
mkdir ".config";
my (@all, @rules);
if(-d ".config"){
    system "rm -r .config";
}
mkdir ".config";
foreach my $dir (@dir_list){
    my @args;
    if($dir=~/\b(\w+)$/){
        my $key;
        if($1 eq "libfabric"){
            $key = "ofi_config_arg";
        }
        elsif($1 eq "openpa"){
            $key = "opa_config_arg";
        }
        else{
            $key = $1."_config_arg";
        }
        if($ENV{$key}){
            push @args, $ENV{$key};
        }
    }
    my @cmds = "cd $dir";
    if(-f "$dir/setup"){
        push @cmds, "source setup";
    }
    my $configure;
    if(-f "$dir/configure.gnu"){
        $configure = "./configure.gnu";
    }
    elsif(-f "$dir/configure"){
        $configure="./configure";
    }
    else{
        die "Missing configure in $dir\n";
    }
    push @cmds, "$configure ". join(' ', @args);
    my $target = $dir;
    $target=~s/\//./g;
    $target=".config/$target";
    push @all, $target;
    push @rules, "$target:";
    push @rules, "\t".join(' && ', @cmds);
    push @rules, "\ttouch $target";
    push @rules, "";
}
open Out, ">.config/Makefile" or die "Can't write .config/Makefile.\n";
print "  --> [.config/Makefile]\n";
my $t = join(' ', @all);
print Out "all: $t\n";
print Out "\n";
foreach my $l (@rules){
    print Out "$l\n";
}
close Out;
system "make -f .config/Makefile -j$opts{j} all" || exit 1;
