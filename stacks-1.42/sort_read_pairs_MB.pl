#!/usr/bin/perl
#
# Copyright 2011-2013, Julian Catchen <jcatchen@uoregon.edu>
#
# This file is part of Stacks.
#
# Stacks is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Stacks is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Stacks.  If not, see <http://www.gnu.org/licenses/>.
#

#
# Sort paired-end sequences according to the stacks the non-paired-end 
# was found in. 
#
# By Julian Catchen <jcatchen@uoregon.edu>
#

use strict;
use constant stacks_version => "_VERSION_";
use Data::Dumper qw(Dumper);
use constant true  => 1;
use constant false => 0;

my $debug          = 0;
my $white_list     = "";
my $cat_white_list = "";
my $cat_white_name = "";
my $in_path        = "";
my $out_path       = "";
my $samp_path      = "";
my $samp_type      = "fastq";
my $out_type       = "fasta";
my $out_read_type       = "single";
my $gzipped        = false;

parse_command_line();

my (@files, %matches, %stacks, %reads, %marker_wl);
build_file_list(\@files);
#~ print "build file_list : ",join(" - ",@files ),"\n";

my ($file, $num_files, $i, $key);
# sauve les id des locus à analyser
if (length($cat_white_list) > 0) {
    load_white_list($cat_white_list, \%marker_wl);
    print STDERR "Loaded ", scalar(keys %marker_wl), " catalog IDs from '$cat_white_list'\n";
    #~ print "whitelist : ",join(" - ", keys %marker_wl) ,":", join(" - ", values %marker_wl),"\n";
}else{
	$marker_wl{$cat_white_name}++;
	print STDERR "Loaded catalog IDs : '$cat_white_name'\n";
	#~ print "whitelist : ",join(" - ", keys %marker_wl) ,":", join(" - ", values %marker_wl),"\n";
}

$num_files = scalar(@files);
$i         = 1;
foreach $file (@files) {
    printf(STDERR "Loading catalog matches, file % 2s of % 2s [%s]\n", $i, $num_files, $file->{'prefix'});
    #
    # Load the sstacks file, listing matches between the stacks and the catalog
    #
    load_matches($in_path, $file, \%matches, \%marker_wl);
	#~ print "matches : ",join(" - ", keys %matches),join(" - ", values %matches), "\n";
    $i++;
}

print STDERR "\n";
#
# Check if files already exist, if so, exit to prevent adding data to existing files.
#
my ($cat_id, $path);

#~ print STDERR join ( " - " , keys %{$matches{"9320"}}) , "\n";
#~ print STDERR join ( " - " , values %{$matches{"9320"}}) , "\n";


foreach $cat_id (keys %matches) {
    #
    # Check that this catalog ID only has a single match from each sample.
    #
    
    if (check_mult_catalog_matches($matches{$cat_id}) == true){
		print STDERR $cat_id." not single match for each sample\n" if $debug;
		#~ next;
	}
	
	# verification que les fichiers de sorties n'existes pas déjà.
    $path  = $out_path . "/" . $cat_id."_2" ;
    $path .= $out_type eq "fasta" ? ".fa" : ".fq";
    if (-e $path) {
	die("Error: output files already exist. This program will append data to files if\n" .
	    "they already exist. Please delete these files and re-execute sort_read_pairs.pl.\n");
    }
}

#~ print STDERR join ( " - " , keys %{$matches{"9320"}}) , "\n";
#~ print STDERR join ( " - " , values %{$matches{"9320"}}) , "\n";

$i = 1;
foreach $file (@files) {
    printf(STDERR "Processing file % 2s of % 2s [%s]\n", $i, $num_files, $file->{'prefix'});
	$i++;
    #
    # Load the ustacks tag file for each sample, containing the read-to-stack mappings
    #
    $stacks{$file->{'prefix'}} = {};

    my @matches_array;
    foreach $cat_id (keys %matches){
		@matches_array=(@matches_array, grep $_ =~ qr/$file->{'prefix'}/i , keys %{$matches{$cat_id}});
	}
	
    if (scalar(@matches_array)==0){
		print STDERR "  no match\n";
		next;
	}
    
    print STDERR "  Loading tag file...\n";
    load_stacks($in_path, $file, $stacks{$file->{'prefix'}}, \@matches_array);
    print STDERR "done.\n";

    #
    # Map the read-pairs to the stack/catalog match they correspond to.
    #
    if ($out_read_type eq "pair"){
		$reads{$file->{'prefix'}} = {};

		print STDERR "  Loading sample file 1... ";
		process_read_pairs($samp_path, $file, \%stacks, $reads{$file->{'prefix'}}, 1);
		print STDERR "done.\n";
		

		print STDERR "  Printing results... ";
		print_results($out_path, \%matches, \%stacks, \%reads, 1);
		print STDERR "done.\n";

		#
		# Clean up memory usage.
		#
		undef(%{$reads{$file->{'prefix'}}});
	}
    
    $reads{$file->{'prefix'}} = {};

    print STDERR "  Loading sample file 2 ... ";
    process_read_pairs($samp_path, $file, \%stacks, $reads{$file->{'prefix'}}, 2);
    print STDERR "done.\n";

    print STDERR "  Printing results... ";
    print_results($out_path, \%matches, \%stacks, \%reads, 2);
    print STDERR "done.\n";

    #
    # Clean up memory usage.
    #
    print STDERR "  Clearing memory...";
    undef(%{$stacks{$file->{'prefix'}}});
    undef(%{$reads{$file->{'prefix'}}});
    print STDERR "  done.\n";
}

sub load_matches {
    my ($in_path, $in_file, $matches, $marker_wl) = @_;

    my ($file, $in_fh, $line, @parts, $key);

    if ($gzipped == true) {
	$file  = $in_path . "/" . $in_file->{'prefix'} . $in_file->{'suffix'} . ".matches.tsv.gz";
	open($in_fh, "gunzip -c $file |") or die("Unable to open catalog matches file '$file', $!\n");
    } else {
	$file  = $in_path . "/" . $in_file->{'prefix'} . $in_file->{'suffix'} . ".matches.tsv";
	open($in_fh, "<$file") or die("Unable to open catalog matches file '$file', $!\n");
    }

    while ($line = <$in_fh>) {
        chomp $line;
        @parts = split(/\t/, $line);

        if (length($cat_white_list) > 0 or $cat_white_name ne "") {
			next if (!defined($marker_wl->{$parts[2]}));
		}

        if (!defined($matches->{$parts[2]})) {
            $matches->{$parts[2]} = {};
        }

        #
        # Index by catalog_ID -> sample_ID|stack_ID
        #
        $key = $in_file->{'prefix'} . "|" . $parts[4];
        $matches->{$parts[2]}->{$key}++;
    }

    close($in_fh);
}

# load_stacks($in_path, $file, $stacks{$file->{'prefix'}}, \@matches_array);
sub load_stacks {
    my ($in_path, $in_file, $stacks, $matches_array) = @_;

    my ($file, $in_fh, $line, @parts);

    if ($gzipped == true) {
	$file = $in_path . "/" . $in_file->{'prefix'} . $in_file->{'suffix'} . ".tags.tsv.gz";
	open($in_fh, "gunzip -c $file |") or die("Unable to open '$file', $!\n");
    } else {
	$file = $in_path . "/" . $in_file->{'prefix'} . $in_file->{'suffix'} . ".tags.tsv";
	open($in_fh, "<$file") or die("Unable to open '$file', $!\n");
    }

    while ($line = <$in_fh>) {
        chomp $line;
        @parts = split(/\t/, $line);

        next if ($parts[6] eq "consensus" || $parts[6] eq "model");
        if (! grep {$_ eq $in_file->{'prefix'}."|".$parts[2]} @{$matches_array}){
			next;	
		}

        #
        # Index by sequence ID -> stack ID
        #
        my $read_id=substr($parts[8], 0, -2);
        my $search= $in_file->{'prefix'}."_1.fq_" ;			# add Maria
        $read_id =~ s/$search//ee ;							# add Maria
        $stacks->{$read_id} = $parts[2];
        #~ print "load_stacks " , $read_id , $parts[2], "\n";
    }

    close($in_fh);
}

# process_read_pairs($samp_path, $file, \%stacks, $reads{$file->{'prefix'}}, 1);
sub process_read_pairs {
    my ($in_path, $in_file, $stacks, $reads, $num_read) = @_;

    my ($file, $in_fh, $line, $seq, $qual, $key, $read_id);
	
    if ($in_file->{'suffix'} eq ".1") {
		if ($samp_type eq "gzfastq") {
			$file = $in_path . "/" . $in_file->{'prefix'} . ".$num_read.fq.gz";
			open($in_fh, "gunzip -c $file |") or die("Unable to open paired-end input file '$file'\n");
		} else {
			$file = $in_path . "/" . $in_file->{'prefix'} . ".$num_read.fq";
			open($in_fh, "<$file") or die("Unable to open paired-end input file '$file'\n");
		}
    } else {
		if ($samp_type eq "gzfastq") {
			$file = $in_path . "/" . $in_file->{'prefix'} . "_$num_read.fq.gz";
			#~ print $file, "\n";
			open($in_fh, "gunzip -c $file |") or die("Unable to open paired-end input file '$file'\n");
		} else {
			$file = $in_path . "/" . $in_file->{'prefix'} . "_$num_read.fq";
			open($in_fh, "<$file") or die("Unable to open paired-end input file '$file'\n");
		}
    }
	
    while ($line = <$in_fh>) {
		next if (substr($line, 0, 1) ne "@");
		chomp $line;

	# version après process_rad_tags
		$read_id = substr($line, 1, -2);
		my $search= $in_file->{'prefix'}."_$num_read.fq_" ;
		$read_id =~ s/$search//ee ;
	# versin sans process_rad_tags
		#~ $line =~ s/ /:/ ;
		#~ my @list = split(':', $line);
		#~ $read_id = join("_",@list[3],@list[4],@list[5],@list[6]) ; 
		#~ print "process_read_pairs " , $read_id , "\n" ;
		$seq     = <$in_fh>;
		chomp $seq;
		
		#
		# Read the repeated ID and the quality scores.
		#
		<$in_fh>;
		$qual = <$in_fh>;
		chomp $qual;

		#~ print "liste des read_id : \n";
		#~ print Dumper $stacks->{$in_file->{'prefix'}} ; 
        $key = $stacks->{$in_file->{'prefix'}}->{$read_id};
        next if (!defined($key));
		#~ print $in_file->{'prefix'} , " : read id ",$read_id," stacks ID ", $key,"\n" ;
        if (!defined($reads->{$key})) {
            $reads->{$key} = [];
        }

        push(@{$reads->{$key}}, {'id' => $read_id, 'seq' => $seq, 'qual' => $qual});
    }
}

sub print_results {
    my ($out_path, $matches, $stacks, $reads, $num_read) = @_;

    my ($path, $cat_id, $sample, $stack_id, $read, $out_fh, $i, @keys, $count, $key, $mult_hits);

    # 
    # If a catalog ID matches stacks from multiple samples, print them out together.
    #
    foreach $cat_id (keys %{$matches}) {
        #~ #
        #~ # Check that this catalog ID only has a single match from each sample.
        #~ #
        #~ next if (check_mult_catalog_matches($matches->{$cat_id}) == true);		# Maria: already done at the begining

        $path  = $out_path . "/" . $cat_id."_".$num_read;
		$path .= $out_type eq "fasta" ? ".fa" : ".fq";
        open($out_fh, ">>$path") or die("Unable to open $path; '$!'\n");

        foreach $key (keys %{$matches->{$cat_id}}) {			# key = sample_ID|stack_ID

            ($sample, $stack_id) = split(/\|/, $key);
			
			foreach $read (@{$reads->{$sample}->{$stack_id}}) {
				if ($out_type eq "fasta") {
					print $out_fh ">", $cat_id, "|", $sample, "|", $stack_id, "|", $read->{'id'}, "\n",
					$read->{'seq'}, "\n";
				} else {
					print $out_fh "@", $cat_id, "|", $sample, "|", $stack_id, "|", $read->{'id'}, "\n",
					$read->{'seq'}, "\n",
					"+\n",
					$read->{'qual'}, "\n";
				}
            }
        }

        close($out_fh);
    }
}

sub check_mult_catalog_matches {
    my ($catalog) = @_;

    my (%samples, @keys, $key, $sample, $stack_id);

	# compte le nombre de stack_id par individus
    foreach $key (keys %{$catalog}) {
        ($sample, $stack_id) = split(/\|/, $key);
        $samples{$sample}++;
    }

	
	# compte le nombre d'échantillon à plusieurs stacks
    my %mult_hits;
    foreach $key (keys %samples) {
        if ($samples{$key} > 1){
			$mult_hits{$key}=1;
		}
    }
    
    #~ print STDERR  "clé multihit: ", join("-", keys %mult_hits),"\n";
    #~ print  STDERR "valeur multihit: ", join("-", values %mult_hits),"\n";
	#~ print STDERR  "clé catalog: ", join("-", keys %{$catalog}),"\n";
	#~ print  STDERR "valeur catalog: ", join("-", values %{$catalog}),"\n"; 
	
    
    # supprime les échantillon à plusieurs match ===> add Maria
    my %test;
    if (keys(%mult_hits) ) {
		foreach $key (keys %{$catalog}) {
			($sample, $stack_id) = split(/\|/, $key);
			if ( defined $mult_hits{$sample} ) {
				delete ($catalog->{$key} );
			} 
		}
        return true;
    } else {
        return false;
    }
}

sub count_reads {
    my ($catalog, $reads) = @_;

    my ($count, $key, $sample, $stack_id);

    $count = 0;

    foreach $key (keys %{$catalog}) {
        ($sample, $stack_id) = split(/\|/, $key);

	if (defined($reads->{$sample}->{$stack_id})) {
	    $count += scalar(@{$reads->{$sample}->{$stack_id}});
	}
    }

    return $count;
}

# appel : build_file_list(\@files); 
# @file est un tableau vide avant appel de build_file_list
sub build_file_list {
    my ($files) = @_;
    my (@ls, $line, $file, $prefix, $suffix);

    # Load a white list of files to process if it is supplied. (les échantillons pas les locus!!)
    my %wl;
    if (length($white_list) > 0) {
		load_white_list($white_list, \%wl);
		print STDERR "Loaded ", scalar(keys %wl), " filenames from '$white_list'\n";
    }

	# in_path correspond au dossier de résultats des fichiers stacks
    @ls = glob("$in_path/*.tags.tsv");
    if (scalar @ls == 0) {
		@ls = glob("$in_path/*.tags.tsv.gz");
		$gzipped = true if (scalar @ls > 0);
    }

    foreach $line (@ls) {
	chomp $line;

	next if (length($line) == 0);	
	next if ($line =~ /batch_\d+\.catalog/);

	($file) = ($line =~ /$in_path\/(.+)\.tags\.tsv\.?g?z?$/); 

    if (length($white_list) > 0) {
	    next if (!defined($wl{$file}));
	}

	
	if ($file =~ /\.1$/) {
		($prefix, $suffix) = ($file =~ /^(.+)(\.1)$/);
	} else {
	    $prefix = $file;
	    $suffix = "";
	}

	#~ print $prefix," ",$suffix, "\n";
	push(@{$files}, {'prefix' => $prefix, 'suffix' => $suffix});
    }
}

sub load_white_list {
    my ($list, $wl) = @_;

    open(WHITE, "<" . $list) 
	or die("Unable to open white list file '$white_list': $!\n");

    my $line   = "";

    while ($line = <WHITE>) {
	chomp $line;

	next if (length($line) == 0);
	next if ($line =~ /^\s*#/);

	$wl->{$line}++;
    }

    close(WHITE);
}

sub parse_command_line {
    while (@ARGV) {
	$_ = shift @ARGV;
	if    ($_ =~ /^-p$/) { $in_path    = shift @ARGV; }
	elsif ($_ =~ /^-o$/) { $out_path   = shift @ARGV; }
	elsif ($_ =~ /^-s$/) { $samp_path  = shift @ARGV; }
	elsif ($_ =~ /^-f$/) { $samp_type  = shift @ARGV; }
	elsif ($_ =~ /^-t$/) { $out_type   = shift @ARGV; }
	elsif ($_ =~ /^-r$/) { $out_read_type   = shift @ARGV; }
	elsif ($_ =~ /^-W$/) { $white_list = shift @ARGV; }
	elsif ($_ =~ /^-w$/) { $cat_white_list = shift @ARGV; }
	elsif ($_ =~ /^-n$/) { $cat_white_name = shift @ARGV; }
	elsif ($_ =~ /^-d$/) { $debug++; }
	elsif ($_ =~ /^-v$/) { version(); exit(); }
	elsif ($_ =~ /^-h$/) { usage(); }
	else {
	    print STDERR "Unknown command line option: '$_'\n";
	    usage();
	}
    }
	if ($cat_white_name ne "" && $cat_white_list ne "") {
	pritn STDERR "give white catalog ID name (-n) OR white catalog ID list (-w).\n";
	usage();
    }

    if ($out_type ne "fasta" && $out_type ne "fastq") {
	pritn STDERR "Output type must be either 'fasta' or 'fastq'.\n";
	usage();
    }
    if ($samp_type ne "fastq" && $samp_type ne "gzfastq") {
	pritn STDERR "Output type must be either 'fastq' or 'gzfastq'.\n";
	usage();
    }
    if ($out_read_type ne "single" && $out_read_type ne "pair") {
	pritn STDERR "Output read type must be either 'single' or 'paire'.\n";
	usage();
    }
    $in_path   = substr($in_path, 0, -1)   if (substr($in_path, -1)   eq "/");
    $out_path  = substr($out_path, 0, -1)  if (substr($out_path, -1)  eq "/");
    $samp_path = substr($samp_path, 0, -1) if (substr($samp_path, -1) eq "/");
}

sub version {
    print STDERR "sort_read_pairs.pl ", stacks_version, "\n";
}

sub usage {
    version();

    print STDERR <<EOQ; 
sort_read_pairs.pl -p path -s path -o path [-t type] [-W white_list] [-w white_list] [-d] [-h]
    p: path to the stacks output files.
    s: path to paired-end sample files.
    f: paired-end sample type, either 'fastq' or 'gzfastq'
    o: path to output the collated FASTA files.
    t: output type, either 'fasta' (default) or 'fastq'.
    r: output read type, either 'single' (default) or 'pair'.
    W: a white list of files to process in the input path.
    w: a white list of catalog IDs to include.
    n: a white name of catalog IDs to include.
    h: display this help message.
    d: turn on debug output.

EOQ

exit(0);
}
