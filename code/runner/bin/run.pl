#!/bin/env perl

use warnings FATAL;
use strict;

use Getopt::Long;
use YAML;
use Cwd qw/getcwd abs_path/;
use File::Basename;

my %DATA   = (em => "/afs/desy.de/user/q/quintera/public/Jets/Lists/Sample_040506e.list",
				    	ep => "/afs/desy.de/user/q/quintera/public/Jets/Lists/Sample_0607p.list");
#my %mc     = (em => "/afs/desy.de/user/q/quintera/public/jets/lists/sample_ari_incl_nc_dis_lowq2_040506e.list",
							#ep => "/afs/desy.de/user/q/quintera/public/jets/lists/sample_ari_incl_nc_dis_lowq2_0607p.list");
my %MC     = (em => "/afs/desy.de/user/q/quintera/public/Jets/Lists/Sample_Ariadne_Low_Q2_NC_DIS_05e.list",
			  			ep => "/afs/desy.de/user/q/quintera/public/Jets/Lists/Sample_Ariadne_Low_Q2_NC_DIS_0607p.list");
my %HERWIG = (em => "/afs/desy.de/user/q/quintera/public/Jets/Lists/Sample_Herwig_PHP_QCD_040506e.list",
							ep => "/afs/desy.de/user/q/quintera/public/Jets/Lists/Sample_Herwig_PHP_QCD_0607p.list");

# Returns data-list to process.
sub parseargs {
	my $opt_em; # electron
	my $opt_ep; # positron
	my $opt_mc;
	my $opt_herwig;
	my $opt_data;
	my $opt_help;
	my $opt_prefix;
	my $opt_macro;
	my $opt_all;
	my $opt_list;

	GetOptions(
			e => \$opt_em,
			p => \$opt_ep,
			mc => \$opt_mc,
			herwig => \$opt_herwig,
			data => \$opt_data,
			all => \$opt_all,
			"oprefix=s" => \$opt_prefix, "o=s" => \$opt_prefix,
			"macro=s" => \$opt_macro,
			"list=s" => \$opt_list,
			help => \$opt_help, h => \$opt_help,
	) or die "$!\n";

	if ($opt_help) {
		print("use: $0 <-e|-p> <--mc|--data> [--macro <ROOT-macro.C>]\n");
		exit 0;
	}

	if ($opt_all) {
		return ("all", $opt_prefix, $opt_macro);
	}


	if (!$opt_list && (!($opt_em || $opt_ep) || !($opt_mc || $opt_data || $opt_herwig))) {
		die "see help for usage\n";
	}

	my $list;
	if ($opt_list) {
		$list = abs_path($opt_list);
	} else {
		if ($opt_mc) {
			$list = $opt_em ? $MC{em} : $MC{ep};
		} elsif ($opt_herwig) {
			$list = $opt_em ? $HERWIG{em} : $HERWIG{ep};
		} else {
			$list = $opt_em ? $DATA{em} : $DATA{ep};
		}
	}

	my $macro = abs_path($opt_macro or "MakeHists.C");

	print "\x1b[38;5;202;1mlist:\x1b[0m $list\n";
	print "\x1b[38;5;202;1mmacro:\x1b[0m $macro\n"; 
	print "\x1b[38;5;202;1mprefix:\x1b[0m $opt_prefix\n"; 

	return ($list, $opt_prefix, $macro);
}


sub main {
	my ($listpath, $prefix, $macro) = @_;

	$listpath =~ /([^\/]+)[.]list/;
	my $outdir = getcwd() . "/${prefix}$1";

	my $listlen = (split(/ /, `wc -l $listpath`))[0];

	print "listpath: $listpath\n";
	print "listlen: $listlen\n";
	print "outdir: $outdir\n";

	my %cut;
	if ($ENV{CUT_JET_MINET}) {
		$cut{jet_minEt} = "export CUT_JET_MINET=$ENV{CUT_JET_MINET}";
	}
	if ($ENV{CUT_LEPTON_MINE}) {
		$cut{lepton_minE} = "export CUT_LEPTON_MINE=$ENV{CUT_LEPTON_MINE}";
	}
	if ($ENV{CUT_LEPTON_MINY}) {
		$cut{lepton_miny} = "export CUT_LEPTON_MINY=$ENV{CUT_LEPTON_MINY}";
	}
	if ($ENV{CUT_LEPTON_MAXY}) {
		$cut{lepton_maxy} = "export CUT_LEPTON_MAXY=$ENV{CUT_LEPTON_MAXY}";
	}

	my $batch = <<EOM;
Universe     = vanilla
Initialdir   = $outdir/
Output       = $outdir/condor/condor.\$(Process).out
Error        = $outdir/condor/condor.\$(Process).err
Log          = $outdir/condor/condor.\$(Process).log
Executable   = /bin/sh
Arguments    = $outdir/runInc.sh \$(Process)
Requirements = OpSysAndVer == "SL6"
Notify_user	 = ivan.pidhurskyi\@desy.de
Queue $listlen
EOM

  my $local_macro = $outdir . '/' . basename($macro);
	my $job = <<EOM;
#!/bin/bash

$cut{jet_minEt}
$cut{lepton_minE}
$cut{lepton_miny}
$cut{lepton_maxy}

export REWEIGHT=$ENV{REWEIGHT}
export REWEIGHT_DATA=$ENV{REWEIGHT_DATA}
export REWEIGHT_RECO=$ENV{REWEIGHT_RECO}

echo `date`
counter=0

for infile in `cat $listpath` ; do
	if [ "\$counter" = "\$1" ]; then  # Run parallel each pt bin
		root.exe $local_macro\\(\\"\$infile\\",\\"$outdir/out\\"\\)
		let counter=counter+1
		break
	fi
	let counter=counter+1
done

echo `date`
EOM

	mkdir $outdir or die $!;
	mkdir "$outdir/out" or die $!;
	mkdir "$outdir/condor" or die $!;
	open(BATCH, '>', "$outdir/batch.job") or die $!;
	open(JOB, '>', "$outdir/runInc.sh") or die $!;
	print BATCH $batch;
	print JOB $job;
	close BATCH;
	close JOB;
	system("cp $macro JetOrange2018.h $outdir") == 0 or die $!;

	my $cwd = getcwd();
	chdir $outdir or die $!;
	system("condor_submit batch.job") == 0 or die $!;
	chdir $cwd or die $!;
}

my ($listpath, $PREFIX, $MACRO) = parseargs();
$MACRO ||= "MakeHist_v5_1.C";
if ($PREFIX) {
	mkdir $PREFIX or die $! unless -d $PREFIX;
	$PREFIX = "$PREFIX/";
}

if ($listpath eq 'all') {
	print "> \e[38;5;3msubmitting em data\e[0m\n";
	main($DATA{em}, $PREFIX, $MACRO);
	print "\n> \e[38;5;3msubmitting ep data\e[0m\n";
	main($DATA{ep}, $PREFIX, $MACRO);
	print "\n> \e[38;5;3msubmitting em MC\e[0m\n";
	main($MC{em}, $PREFIX, $MACRO);
	print "\n> \e[38;5;3msubmitting ep MC\e[0m\n";
	main($MC{ep}, $PREFIX, $MACRO);
	print "\n> \e[38;5;3msubmitting em HERWIG\e[0m\n";
	main($HERWIG{em}, $PREFIX, $MACRO);
	print "\n> \e[38;5;3msubmitting ep HERWIG\e[0m\n";
	main($HERWIG{ep}, $PREFIX, $MACRO);
} else {
	main($listpath, $PREFIX, $MACRO);
}

0
