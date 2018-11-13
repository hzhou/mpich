#!/usr/bin/perl
use strict;
our %subst;
our @all;
our @rules;
our %dependency;

sub push_rule_target_dependency {
    my ($t) = @_;
    my $dep = $dependency{"_"};
    $dep .= $dependency{$t};
    $dep=~s/^\s+//;
    if($dep){
        push @rules, "$t: $dep";
    }
    else{
        push @rules, "$t:";
    }
}

my %opts=("-j"=>2);
$opts{skip_binding}=1;
foreach my $a (@ARGV){
    if($a =~ /-j(\d+)/){
        $opts{"-j"}=$1;
    }
}
push @rules, "maint/Version: maint/Version.base.m4 maint/version.m4";
push @rules, "\tautom4te -l M4sugar maint/Version.base.m4 >maint/Version";
push @rules, "";
push @all, "maint/Version";
my $t = `pwd`;
chomp $t;
$subst{abs_mpichsrcdir} = $t;
$subst{abs_srcdir} = "$t/maint";
my $t = `which perl`;
chomp $t;
$subst{PERL} = $t;
my $t = `echo ""| xargs ls | wc -l`;
$t=~s/\s+//g;
if($t ne '0'){
    $subst{XARGS_NODATA_OPT}="-r";
}
else{
    $subst{XARGS_NODATA_OPT}="";
}
$subst{configure_input}="";
my @files=qw(checkbuilds getcoverage genstates clmake extractstrings extractcvars extractstates extractfixme createcoverage gcovmerge createhtmlindex getfuncstack);
push @files, "cvardirs";
foreach my $f (@files){
    my @lines;
    open In, "maint/$f.in" or die "Can't open maint/$f.in.\n";
    while(<In>){
        if(/\@\w+\@/){
            while(/\@(\w+)\@/g){
                if(!defined $subst{$1}){
                    warn "$f: subst: $1 not defined\n";
                }
            }
            s/\@(\w+)\@/$subst{$1}/ge;
        }
        push @lines, $_;
    }
    close In;
    open Out, ">maint/$f" or die "Can't write maint/$f.\n";
    foreach my $l (@lines){
        print Out $l;
    }
    close Out;
    system "chmod a+x maint/$f";
}
open Out, ">src/mpi/errhan/defmsg.h" or die "Can't write src/mpi/errhan/defmsg.h.\n";
print Out "typedef struct { const unsigned int sentinal1; const char *short_name, *long_name; const unsigned int sentinal2; } msgpair;\n";
print Out "static const int generic_msgs_len = 0;\n";
print Out "static msgpair generic_err_msgs[] = { {0xacebad03, 0, \"no error catalog\", 0xcb0bfa11}, };\n";
print Out "static const int specific_msgs_len = 0;\n";
print Out "static msgpair specific_err_msgs[] = {  {0xacebad03,0,0,0xcb0bfa11}, };\n";
print Out "#if MPICH_ERROR_MSG_LEVEL > MPICH_ERROR_MSG__NONE\n";
print Out "#define MPIR_MAX_ERROR_CLASS_INDEX 54\n";
my $t = '0, 'x53 . '0';
print Out "static int class_to_index[] = { $t };\n";
print Out "#endif\n";
close Out;
my $t = ".autogen/extractstates";
push_rule_target_dependency($t);
push @rules, "\t\@echo Creating the enumeration of logging states into src/include/mpiallstates.h";
push @rules, "\tperl maint/extractstates";
push @all, $t;
push @rules, "\ttouch $t";
push @rules, "";
my @tlist = qw(src/mpi src/mpi_t src/nameserv src/util src/binding src/include src/mpid src/pmi src/mutex);
my $t = ".autogen/extractcvars";
push_rule_target_dependency($t);
push @rules, "\t\@echo Extracting control variables";
my $cvardirs=join(' ', @tlist);
push @rules, "\tperl maint/extractcvars --dirs=\"$cvardirs\"";
push @all, $t;
push @rules, "\ttouch $t";
push @rules, "";
$dependency{"subsys_include.m4"}="maint/gen_subcfg_m4";
$dependency{".autogen/main"}.=" subsys_include.m4";
push_rule_target_dependency("subsys_include.m4");
push @rules, "\tperl maint/gen_subcfg_m4";
push @rules, "";
my @tlist=(
    "src/mpl",
    "src/pm/hydra",
    "src/mpi/romio",
    );
foreach my $a (@tlist){
    if(-l "$a/confdb"){
        system "rm $a/confdb";
    }
    my $t2 = $a;
    $t2=~s/\//./g;
    my $t = "$t2.confdb";
    $dependency{".autogen/$t2"}.=" .autogen/$t";
    my $t = ".autogen/$t";
    push_rule_target_dependency($t);
    push @rules, "\tif test -d $a/confdb; then rm -r $a/confdb; fi";
    push @rules, "\tcp -pPR confdb $a/confdb";
    push @rules, "\tcp -pPR maint/version.m4 $a/";
    if($a ne "src/mpl"){
        push @rules, "\tcp -pPR maint/version.m4 $a/version.m4";
    }
    push @rules, "\ttouch $t";
    push @rules, "";
}
my @all_autogens=(
    "src/openpa",
    "src/hwloc",
    "src/izem",
    "src/mpid/ch4/netmod/ucx/ucx",
    "src/mpid/ch4/netmod/ofi/libfabric",
);
foreach my $a (@all_autogens){
    my $t=$a;
    $t=~s/\//./g;
    my $t = ".autogen/$t";
    push_rule_target_dependency($t);
    push @all, $t;
    push @rules, "\t\@echo Running autogen.sh in $a";
    push @rules, "\tcd $a && ./autogen.sh";
    push @rules, "\ttouch $t";
    push @rules, "";
}
my $t = ".autogen/main";
push_rule_target_dependency($t);
push @all, $t;
push @rules, "\t\@echo Running autoreconf in .";
push @rules, "\tautoreconf -ifv";
push @rules, "\ttouch $t";
push @rules, "";
my @all_other=(
    "src/mpl",
    "src/pm/hydra",
    "src/mpi/romio",
    "src/util/logging/rlog",
);
my $t = ".autogen/romio_gluecode";
push_rule_target_dependency($t);
push @rules, "\tcd src/glue/romio && perl all_romio_symbols ../../mpi/romio/include/mpio.h.in";
$dependency{".autogen/src.mpi.romio"}.=" $t";
push @rules, "\ttouch $t";
push @rules, "";
foreach my $a (@all_other){
    my $t=$a;
    $t=~s/\//./g;
    my $t = ".autogen/$t";
    push_rule_target_dependency($t);
    push @all, $t;
    push @rules, "\t\@echo Running autoreconf in $a";
    push @rules, "\tcd $a && autoreconf -ifv";
    push @rules, "\ttouch $t";
    push @rules, "";
}
system "touch src/binding/fortran/mpif_h/Makefile.mk";
system "touch src/binding/fortran/use_mpi/Makefile.mk";
system "touch src/binding/fortran/use_mpi_f08/Makefile.mk";
system "touch src/binding/cxx/mpicxx.h.in";
system "touch src/binding/fortran/mpif_h/mpif.h.in";
system "touch src/binding/fortran/use_mpi/mpi_base.f90.in";
system "touch src/binding/fortran/use_mpi/mpi_constants.f90.in";
system "touch src/binding/fortran/use_mpi_f08/mpi_f08_compile_constants.f90.in";
system "touch src/binding/fortran/use_mpi_f08/mpi_c_interface_types.f90.in";
if(-d ".autogen"){
    system "rm -r .autogen";
}
system "mkdir .autogen";
open Out, ">.autogen/Makefile" or die "Can't write .autogen/Makefile.\n";
print "  --> [.autogen/Makefile]\n";
my $t = join(' ', @all);
print Out "all: $t\n";
print Out "\n";
foreach my $l (@rules){
    print Out "$l\n";
}
close Out;
my $j = "-j".$opts{"-j"};
system "make -f .autogen/Makefile $j all" || exit 1;
