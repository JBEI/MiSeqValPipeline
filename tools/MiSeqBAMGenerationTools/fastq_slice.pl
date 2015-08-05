#!/usr/bin/env perl
use Modern::Perl;
use Getopt::Long;

# DEPENDS ON PIGZ

# divides input fastq into -number reads per file, changing suffixes to @xx.fq.gz for each batch
# also separates read1 / read2 for paired ( -split ) changing suffix to .[12]@xx.fq.gz or .[12].fq.gz
# depending on if sliced into -number batches
#
# TODO: switch file opens to zopen in GT::Utils;
#
my($filename, $reads, $suffix, $help, $prefix, $threads, $do_gzip, $do_split, $ill2std, $verbose);

$suffix = 'fq';
$reads = 2000000;
$threads = 2;

GetOptions(
	'file=s' => \$filename,
    'parts|number=i' => \$reads,
    'gzip!' => \$do_gzip,
    'split!' => \$do_split, #separate read1 read2
    'threads=i' => \$threads,
    'suffix=s' => \$suffix,
    'ill2std|ill2sanger' => \$ill2std, 
    'prefix=s' => \$prefix, #prepended to output filename(s)
    'verbose' => \$verbose,
    'help' => \$help,
) or die $!;


if (!defined $filename || ($filename ne '-' && !-e $filename) || !defined $prefix || $help) {
    print STDERR "\n$0 options input.fq output-prefix\n\n";
    print STDERR " -num [$reads] # reads for each part\n";#
    print STDERR " -prefix [] name before the \@,\n";
    print STDERR " -suffix [$suffix] suffix for filename created\n\n";
    print STDERR " -gzip [N/A] gzip the output files\n";
    print STDERR " -split [N/A] split into read1 read2 ( prefix.1 prefix.2 )\n";
    print STDERR " -threads 2 number of threads for the gzip, it gets io bound especially if reading from a .gz\n";
    print STDERR " more isn't always better\n";
    exit;
}

my $fh;
if ($filename =~ /bz2$/) {
    open($fh, "bzip2 -cd $filename |") or die $!;
}
elsif($filename =~ /gz$/) {
    open($fh, "gzip -cd $filename |") or die $!;
}
elsif($filename eq '-') {
    $fh = * STDIN;
} else {
    open($fh, $filename) or die $!;
}
my $outname1;# = "$prefix\@1.$suffix";
my $outpipe = $ill2std ? "| fqc ill2std - - " : '';
my($oh1, $oh2);

if ($do_gzip) {
    my $out_put = "$outpipe | pigz -p $threads -c ";
    if ($do_split) {
        my $outname1 = "$prefix.1\@1.$suffix";
        my $outname2 = "$prefix.2\@1.$suffix";
        if ($verbose) {
            print "opening $out_put > $outname1.gz\n";
        }
        if ($verbose) {
            print "opening $out_put > $outname2.gz\n";
        }
        open($oh1, "$out_put > $outname1.gz") or die $!;
        open($oh2, "$out_put > $outname2.gz") or die $!;
    } else {
        my $outname = "$prefix\@1.$suffix";
        if ($verbose) {
            print "opening $out_put > $outname.gz\n";
        }
        open($oh1, "$out_put > $outname.gz") or die $!;
        $oh2 = $oh1;
    }
} else {
    if ($do_split) {
        my $outname1 = "$prefix.1\@1.$suffix";
        my $outname2 = "$prefix.2\@1.$suffix";
        warn $outname1;
        if ($outpipe) {
            if ($verbose) {
                print "opening $outpipe > $outname1\n";
            }
            if ($verbose) {
                print "opening $outpipe > $outname2\n";
            }
            open($oh1, "$outpipe > $outname1") or die $!;
            open($oh2, "$outpipe > $outname2") or die $!;
        } else {
            if ($verbose) {
                print "opening > $outname1\n";
            }
            if ($verbose) {
                print "opening > $outname2\n";
            }
            open($oh1, '>', $outname1) or die $!;
            open($oh2, '>', $outname2) or die $!;
        }
    } else {
        my $outname = "$prefix\@1.$suffix";
        if ($outpipe) {
            if ($verbose) {
                print "opening $outpipe > $outname\n";
            }
            open($oh1, "$outpipe > $outname") or die $!;
        } else {
            if ($verbose) {
                print "opening > $outname\n";
            }
            open($oh1, '>', $outname) or die $!;
        }
        $oh2 = $oh1;
    }
}
sub inc_filehandles {
    my($inc, $roh1, $roh2) = @_;
    my $name1 = $do_split ? "$prefix.1\@$inc.$suffix" : "$prefix\@$inc.$suffix";
    my $name2 = "$prefix.2\@$inc.$suffix";
    if ($do_gzip) {
        my $out_put = "$outpipe | pigz -p $threads -c ";
        open($$roh1, "$out_put > $name1.gz") or die "ERROR: can't open $name1.gz for writing from pigz$!";
        if ($do_split) {
            open($$roh2, "$out_put > $name2.gz") or die "ERROR: can't open $name2.gz for writing from pigz$!";
        }
    } else {
        if ($outpipe) {
            open($$roh1, "$outpipe > $name1") or die "ERROR: can't open $name1 for writing $!";
            if ($do_split) {
                open($$roh2, "$outpipe > $name2") or die "ERROR: can't open $name2 for writing $!";
            }
        } else {
            open($$roh1, '>', $name1) or die "ERROR: can't open $name1 for writing $!";
            if ($do_split) {
                open($$roh2, '>', $name2) or die "ERROR: can't open $name2 for writing $!";
            }
        }
    }
}
my $i = 0;
my $output = '';
while ( <$fh> ) {
    if ($i && !($i / 4 % $reads)) {
        inc_filehandles($i / 4 + 1, \$oh1, \$oh2);
    }
    $output.= $_;
    for (0..2) {
        $output.= <$fh> ;
    }
    print $oh1 $output;
    # print STDERR "oh1\n",$output;
    $output = '';
    if ($do_split) {
        for (0..3) {
            $output.= <$fh> ;
        }
        #print STDERR "oh2\n",$output;
        print $oh2 $output;
        $output = '';
    }
    $i += 4;
}
close($fh) or die "can't close fh after reading $!";
close($oh1) or die "can't close fh after writing $!";
if ($do_split) {
    close($oh2) or die "can't close 2nd fh after writing $!";
}