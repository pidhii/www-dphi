#!/bin/env perl

use warnings FATAL;
use strict;

use Getopt::Long;
use Term::ReadKey;
use File::Basename;

my $opt_help;
my $opt_oprefix;
my $opt_verbose;
GetOptions(
		help => \$opt_help, h => \$opt_help,
		"out=s" => \$opt_oprefix, "o=s" => \$opt_oprefix,
    verbose => \$opt_verbose, v => \$opt_verbose,
) or die "$!\n";

if ($opt_help) {
	print "use: $0 [--out <prefix>] <data>\n";
	exit 0;
}

die "see `$0 --help` for usage\n" unless $ARGV[0];

my @DIRS = @ARGV;
my $OPREFIX = $opt_oprefix ? "$opt_oprefix/" : "";

print map { "> \"$_/out/*.root\" -> \"${OPREFIX}" . (fileparse($_))[0] . ".root\"\n" } @DIRS;
print "> confirm [yn]? ";
while (1) {
	ReadMode 'cbreak';
	my $ch = getc(STDIN);
	ReadMode 'normal';
	print "$ch\n" and last if $ch =~ /y/i;
	print "$ch\n" and exit(0) if $ch =~ /n/i;
}

mkdir $OPREFIX or die $! unless -d $OPREFIX;

my @pids;
$| = 1;
foreach my $dir (@DIRS) {
	print "> starting extraction of \"$dir\"\n";
	my $outname = fileparse($dir);
	$outname = "$OPREFIX$outname.root";
	my $pid = fork();
	if ($pid == 0) {
    if ($opt_verbose) {
      system("hadd -k $outname $dir/out/*.root &>/dev/null") == 0 or die $!;
    } else {
      system("hadd -k $outname $dir/out/*.root") == 0 or die $!;
    }
		print "> \"$outname\" is \e[38;5;2mready\e[0m\n";
		exit 0;
	} else {
		push @pids, $pid;
	}
}
$| = 0;
waitpid($_, 0) foreach @pids;

0
