#!/usr/bin/env perl

use Getopt::Long;
use Pod::Usage;
use XML::Simple;
use Spreadsheet::WriteExcel;
use English qw(-no_match_vars);
use Carp;
use CGI qw(:standard);
use POSIX qw(strftime);
use Bio::GeneDesign;

use strict;
use warnings;

my $VERSION = '1.00';
my $JGISBV = "JGISB_Analyze_PacBio_$VERSION";

local $OUTPUT_AUTOFLUSH = 1;

##Get Arguments
my %p = ();
GetOptions (
      'help'          => \$p{HELP},
      'analysisids=s' => \$p{INPUT},
);
pod2usage(-verbose=>99) if ($p{HELP});

################################################################################
################################ SANITY  CHECK #################################
################################################################################
my $GD = Bio::GeneDesign->new();
my $linkbase = 'http://www.broadinstitute.org/igv/projects/current/igv.php?';
$linkbase .= 'sessionURL=http://synbio.jgi-psf.org/pacbio/analyses/igvs/';
#/pacbioSB
my $xllink = 'http://synbio.jgi-psf.org/pacbio/analyses/excels/';
my $base = '/projectb/projectdirs/RD/synbio/pacbioSB/testing/';
#/projectb/projectdirs/RD/synbio/pacbioSB/
my $pathbase = $base . 'jobs/';
my $refbase = $base . 'references/';

#The input file must exist and be a format we care to read.
die "\n JGISB_ERROR: You must supply an analysis ID string.\n"
  if (! $p{INPUT});

my %RANKS = (
  'pass' => 0,
  'fixable' => 1,
  'poorfix' => 2,
  'fail: unfixable' => 6,
  'fail: varcount' => 7,
  'fail: depth' => 8,
  'fail: coverage' => 9,
);
################################################################################
################################# CONFIGURING ##################################
################################################################################
my %ALL;
my %REFS;
my %xREFS;
my %REFLENS;
my %REFOBJS;

my @ALIST = split(q{,}, $p{INPUT});

my $coveragecount = 0;
my $failcount = 0;
my $fixcount = 0;
my $passcount = 0;

################################################################################
################################## ANALYZING ###################################
################################################################################
foreach my $aid (@ALIST)
{
  my $apath = $base . $aid . q{/};
  my $xlpath = $apath . $aid . '.xls';
  my $arefpath = $apath . 'splits/';
  my $xmlpath = $apath . 'igvs/';
  if (! -e $arefpath)
  {
    mkdir $arefpath;
    system "chmod -R 0777 $arefpath";
  }
  if (! -e $xmlpath)
  {
    mkdir $xmlpath;
    system "chmod -R 0777 $xmlpath";
  }

  # GATHER REFERENCES
  #
  #
  my $refpath = $apath . 'references.fasta';
  my ($iter, $filename, $suffix) = $GD->import_seqs($refpath);
  while ( my $obj = $iter->next_seq() )
  {
    my $pid = $obj->id;
    #my $preref = $obj->id;
    #my ($pref, $pid) = split m{\|}x, $preref;
    my $seq = $obj->seq;
    $obj->id($pid);
    $REFS{$pid} = $seq;
    #$REFSx{$pref} = $pid;
    #$xREFS{$pid} = $pref;
    my $reflen = length $seq;
    $REFLENS{$pid} = $reflen;
    #$REFOBJS{$pid} = [$pref, $obj];
    $REFOBJS{$pid} = $obj;
  }
  my @REFLIST = sort keys %REFS;
  my %DATA = ();

  # POOLS
  #
  #
  opendir(my $adir, $apath) || die "Can't opendir $apath: $OS_ERROR\n";
  my @list = grep {$_ !~ m{^\.\.?$}} readdir($adir);
  closedir $adir;
  my @plist = grep {$_ =~ m{pool\d}} @list;
  foreach my $pool (@plist)
  {
    print "Working on $pool\n";
    $DATA{$pool} = {};

    ################################ GATHERING #################################
    my $ppath = $apath . $pool . q{/};
    opendir (my $pdir, $ppath) || die "Can't opendir $ppath: $OS_ERROR\n";
    my @alnlist = grep {$_ =~ m{^\d+$}} readdir($pdir);
    closedir $pdir;
    foreach my $aln (@alnlist)
    {
      my $alnpath = $ppath . $aln . q{/};
      print "\tLooking at alignment $aln\n";
      my $ccs = 0;
      my $settingsfile = $alnpath . 'settings.xml';
      my @settingsarr = slurp($settingsfile);
      my @flags = grep {$_ =~ m{CCS}x} @settingsarr;
      $ccs = 1 if (scalar @flags);

      if ($ccs == 1)
      {
        my $gatkpath = $alnpath . 'GATK/';

        # Gather GATK coverage
        #
        #
        my $gcovpath = $gatkpath . 'covdepth';
        my @gcovlines = slurp($gcovpath);
        shift @gcovlines;
        my %gatkseens = map {$_ => []} keys %REFS;
        my %gatkdepth = map {$_ => 0} keys %REFS;
        foreach my $line (@gcovlines)
        {
          my ($namebase, $cov, $avd, $cont) = split m{\s}, $line;
          my ($name, $base) = split q{:}, $namebase;
          push @{$gatkseens{$name}}, $cov if ($cov > 0);
          $gatkdepth{$name} += $cov;
        }
        foreach my $pid (keys %REFS)
        {
          my $depth = int ($gatkdepth{$pid} / $REFLENS{$pid});
          $DATA{$pool}->{$pid}->{GATK}->{depth} = $depth;
          my @lines = @{$gatkseens{$pid}};
          if ($REFLENS{$pid} == scalar @lines)
          {
            $DATA{$pool}->{$pid}->{GATK}->{coverage} = 100.00;
          }
          else
          {
            $DATA{$pool}->{$pid}->{GATK}->{coverage} = 0;
          }
          $DATA{$pool}->{$pid}->{GATK}->{variance} = [];
        }

        # Gather GATK variance
        #
        #
        my $gvarpath = $gatkpath . 'snps.gatk.vcf';
        my @gvlines = slurp($gvarpath);
        foreach my $line (@gvlines)
        {
          next if $line =~ m{\#}x;
          my @atts = split m{\s}x, $line;
          my $pid = $atts[0];
          if (! exists $REFS{$pid})
          {
            print "\t\t$pid not in reference!\n";
          }
          push @{$DATA{$pool}->{$pid}->{GATK}->{variance}}, $line;
        }
      }

      elsif ($ccs == 0)
      {
        # Gather GCON coverage
        #
        #
        my $ocovpath = $alnpath . 'results/variants.xml';
        my $xml = new XML::Simple;
        my $data = $xml->XMLin($ocovpath);
        my $atts = $data->{attributes}->{attribute};
        my %gconseens;
        while ( my ($key, $value) = each(%{$atts}) )
        {
          my $val = $value->{value} || 0.00;
          my ($leg, $pid) = split q{ - }, $key;
          $gconseens{$pid}++;
          if (! exists $REFS{$pid})
          {
            print "\t\t$pid not in reference!\n";
          }
          if ($leg eq 'Bases Called (%)')
          {
            $DATA{$pool}->{$pid}->{GCON}->{coverage} = $val;
            $DATA{$pool}->{$pid}->{GCON}->{variance} = [];
          }
        }

        # Gather GCON variance
        #
        #
        my $ovarpath = $alnpath . 'data/variants.vcf';
        my @ovlines = slurp($ovarpath);
        foreach my $line (@ovlines)
        {
          next if $line =~ m{\#}x;
          my @atts = split m{\s}x, $line;
          my $pid = $atts[0];
          if (! exists $REFS{$pid})
          {
            print "\t\t$pid not in reference!\n";
          }
          push @{$DATA{$pool}->{$pid}->{GCON}->{variance}}, $line;
          #print $line, "\n";
        }
      }
    }

    ################################# GRADING ##################################
    my %STATUS = ();
    my %LINKS = ();
    my %LOCI = ();
    my %FIXES = ();
    my @PIDS = keys %{$DATA{$pool}};
    foreach my $pid (@PIDS)
    {
      #print "\tEvaluating $pid:\n";
      my $runpath = $xmlpath . $pid . ".xml";
      my $runbase = $linkbase . $pool . q{/} . $pid . ".xml";
      open (my $IGVXML, '>', $runpath) || croak "Can't open $runpath: $OS_ERROR\n";
      print $IGVXML xmly($pool, $pid);
      close $IGVXML;
      $LINKS{$pid} = $runbase;

      # Coverage
      #
      #
      my $coverage = 0.00;
      $coverage += $DATA{$pool}->{$pid}->{GCON}->{coverage};
      $coverage += $DATA{$pool}->{$pid}->{GATK}->{coverage};
      if ($coverage < 200.00)
      {
        $STATUS{$pid} = 'fail: coverage';
        #print "\t\tfailing (coverage) $pid\n";
        $coveragecount++;
        next;
      }

      # Coverage depth
      #
      #
      #my $depth = $DATA{$pool}->{$pid}->{GATK}->{depth};
      #if ($depth < 10)
      #{
      #  $STATUS{$pid} = 'fail: depth';
      #  #print "\t\tfailing (depth of coverage) $pid\n";
      #  $coveragecount++;
      #  next;
      #}


      # Variants
      #
      #
      my @gconvars = @{$DATA{$pool}->{$pid}->{GCON}->{variance}};
      my @gatkvars = @{$DATA{$pool}->{$pid}->{GATK}->{variance}};
      my $gconcount = scalar @gconvars;
      my $gatkcount = scalar @gatkvars;
      my $varcount = $gconcount + $gatkcount;
      if ($varcount == 0)
      {
        $DATA{$pool}->{$pid}->{RESULT} = 'pass';
        $STATUS{$pid} = 'pass';
        #print "\t\tpassing $pid\n";
        $passcount++;
        next;
      }
      elsif ($varcount > 2)
      {
        $STATUS{$pid} = 'fail: varcount';
        #print "\t\tfailing (variant count) $pid\n";
        next;
      }
      elsif ($gconcount == 0)
      {
        @gconvars = @gatkvars;
      }
      elsif ($gatkcount == 0)
      {
        @gatkvars = @gconvars;
      }
      my @fixable = compare_variants(\@gatkvars, \@gconvars);

      if (! scalar @fixable)
      {
        $STATUS{$pid} = 'fail: unfixable';
        #print "\t\tfailing (unfixable) $pid\n";
        $failcount++;
        next;
      }
      my @coords = ($fixable[0], $fixable[1]);
      my ($bit, $upol, $dnol) = fixing_primers(\@coords, $REFS{$pid});
      if (! defined $bit)
      {
        $STATUS{$pid} = 'fail: unfixable';
        #print "\t\tfailing (unfixable) $pid\n";
        $failcount++;
        next;
      }
      if ($bit > 1)
      {
        $STATUS{$pid} = 'poorfix';
        #print "\t\tpoorfixing $pid\n";
        $fixcount++;
      }
      else
      {
        $STATUS{$pid} = 'fixable';
        #print "\t\tfixing $pid\n";
        $fixcount++;
      }
      $LOCI{$pid} = $fixable[2];
      #print "\t\tWill try to fix $pid from (@fixable)\n";
      $FIXES{$pid} = [$upol, $dnol];
    }
    $ALL{$pool}->{STATUS} = \%STATUS;
    $ALL{$pool}->{FIXES} = \%FIXES;
    $ALL{$pool}->{LOCI} = \%LOCI;
    $ALL{$pool}->{LINKS} = \%LINKS;
  }

  foreach my $pid (@REFLIST)
  {
    my $obj = $REFOBJS{$pid};
    my $refname = $pid;
    $GD->export_seqs(
      -filename => $pid . '.fasta',
      -path => $arefpath,
      -format => 'fasta',
      -sequences => [$obj]
    );
    my $faipath = $arefpath . $pid . '.fasta.fai';
    my $offset = 2 + length $refname;
    open (my $FAI, '>', $faipath) || croak "Can't open $faipath: $OS_ERROR\n";
    print $FAI $refname, q{ }, $REFLENS{$pid}, q{ }, $offset, q{ }, 60, q{ }, 61, "\n";
    close $FAI;
  }

  my $totalcount = $fixcount + $failcount + $passcount + $coveragecount;
  $totalcount = $totalcount || 0.1;
  my $statistic = q{};
  my $fixperc = sprintf("%.1f", (100 * $fixcount / $totalcount));
  my $failperc = sprintf("%.1f", (100 * $failcount / $totalcount));
  my $passperc = sprintf("%.1f", (100 * $passcount / $totalcount));
  my $covperc = sprintf("%.1f", (100 * $coveragecount / $totalcount));

  ################################ HTML REPORT #################################
  my %BESTS;
  if (1)
  {
    my $style = <<"END";
<!--
td,th
{
  text-align: center;
}
-->
END

    my %COLORS = (
      'pass' => "\#0F0",
      'fixable' => "\#FF0",
      'poorfix' => "\#FA0",
      'fail: unfixable' => "\#F00",
      'fail: varcount' => "\#B00",
      'fail: depth' => "\#EEE",
      'fail: coverage' => "\#AAA",
    );
    my $html = start_html(-style => {-code => $style});
    my $xllinky = a({href => $xllink}, 'excel file');
    $html .= "<br>$xllinky<br>";
    $html .= "<table border=\"1\" align=\"center\"><tr>";
    $html .= "<td>Passing: $passperc\%</td>";
    $html .= "<td>Fixable: $fixperc\%</td>";
    $html .= "<td>Failing: $failperc\%</td>";
    $html .= "<td>Bad coverage: $covperc\%</td>";
    $html .= "</tr></table>";
    $html .= "<table border=\"1\" align=\"center\">\n";
    $html .= "<tr>\n\t<td>color key</td>\n";
    foreach my $status (sort {$RANKS{$a} <=> $RANKS{$b}} keys %COLORS)
    {
      my $bgcolor = $COLORS{$status};
      $html .= "<td bgcolor=\"$bgcolor\">$status</td>\n";
    }
    $html .= "</tr>\n";
    $html .= "</table>\n";
    $html .= "<table border=\"1\" align=\"center\">\n";
    $html .= "<tr>\n\t<td>clone</td>\n";
    foreach my $jid (@plist)
    {
      $html .= "\t<th>$jid</th>\n";
    }
    $html .= "\t<td> Matt's Best Bet</td>\n";
    $html .= "</tr>\n";
    foreach my $pid (@REFLIST)
    {
      my $best = ['fail: coverage', q{_}, 0];
      $html .= "<tr>\n";
      $html .= "\t<th> $pid</th>\n";
      foreach my $jid (@plist)
      {
        my $stat = $ALL{$jid}->{STATUS}->{$pid};
        my $link = $ALL{$jid}->{LINKS}->{$pid};
        my $depth = $DATA{$jid}->{$pid}->{GATK}->{depth};
        my $bgcolor = $COLORS{$stat};
        my $text = a({href => $link}, $depth);
        my $pname = a({href => $link}, $jid . " ($depth)");
        my $rank = $RANKS{$stat};
        my $beat = $RANKS{$best->[0]};
        if ($rank < 5 && (($rank < $beat && $depth > 10) || ($rank <= $beat && $depth > $best->[2])))
        {
          $BESTS{$pid} = $jid;
          $best->[0] = $stat;
          $best->[1] = $pname;
          $best->[2] = $depth;
        }
        $html .= "\t<td bgcolor=\"$bgcolor\">$text</td>\n";
      }
      my $allcolor = $COLORS{$best->[0]};
      $html .= "\t<td bgcolor=\"$allcolor\">" . $best->[1] . "</td>\n";
      $html .= "</tr>";
    }
    $html .= "</table>";
    $html .= end_html();
    my $htmlpath = $apath . $aid . '.html';
    open (my $REPORT, '>', $htmlpath) || croak "Can't open $htmlpath: $OS_ERROR\n";
    print $REPORT $html;
    close $REPORT;
    system "chmod 0777 $htmlpath";
  }


  ############################### EXCEL  REPORT ################################
  if (1)
  {
    my $track = Spreadsheet::WriteExcel->new($xlpath);
    $track->compatibility_mode();

    my $header = $track->add_format();
    $header->set_bold();
    my $pass = $track->add_format();
    $pass->set_bg_color('green');
    my $gfix = $track->add_format();
    $gfix->set_bg_color('yellow');
    my $pfix = $track->add_format();
    $pfix->set_bg_color(52);
    my $fail = $track->add_format();
    $fail->set_bg_color('red');
    my $ufix = $track->add_format();
    $ufix->set_bg_color(29);
    my $coverage = $track->add_format();
    $coverage->set_bg_color('grey');
    my $covdep = $track->add_format();
    $covdep->set_bg_color('47');

    my %XLCOLORS = (
      'pass' => $pass,
      'fixable' => $gfix,
      'poorfix' => $pfix,
      'fail: coverage' => $coverage,
      'fail: depth' => $covdep,
      'fail: varcount' => $fail,
      'fail: unfixable' => $ufix,
    );


    my $bestsheet = $track->add_worksheet('bestbets');
    $bestsheet->write(0, 0, 'clone', $header);
    $bestsheet->write(0, 1, 'best pool', $header);
    $bestsheet->write(0, 2, 'variant', $header);
    $bestsheet->write(0, 3, 'univF - reverse primer', $header);
    $bestsheet->write(0, 4, 'forward - univR primer', $header);
    my $mastersheet = $track->add_worksheet("overview");
    my $x = 1;
    foreach my $jid (@plist)
    {
      my $currsheet = $track->add_worksheet($jid);
      $mastersheet->write(0, $x, $jid, $header);
      $currsheet->write(0, 0, 'clone', $header);
      $currsheet->write(0, 1, 'result', $header);
      $currsheet->write(0, 2, 'variant', $header);
      $currsheet->write(0, 3, 'univF - reverse primer', $header);
      $currsheet->write(0, 4, 'forward - univR primer', $header);
      my $y = 1;
      foreach my $pid (@REFLIST)
      {
        $currsheet->write($y, 0, $pid);
        if ($x == 1)
        {
          $mastersheet->write($y, 0, $pid);
          $bestsheet->write($y, 0, $pid);
        }
        my $stat = $ALL{$jid}->{STATUS}->{$pid};
        $currsheet->write($y, 1, $stat, $XLCOLORS{$stat});
        $mastersheet->write($y, $x, $jid . q{ : } . $stat, $XLCOLORS{$stat});
        if (exists $BESTS{$pid} && $BESTS{$pid} eq $jid)
        {
          $bestsheet->write($y, 1, $jid, $XLCOLORS{$stat});
        }
        if ($stat eq 'fixable' || $stat eq 'poorfix')
        {
          my $loci = $ALL{$jid}->{LOCI}->{$pid};
          $currsheet->write($y, 2, $loci, $XLCOLORS{$stat});
          my $fixes = $ALL{$jid}->{FIXES}->{$pid};
          $currsheet->write($y, 3, $fixes->[0], $XLCOLORS{$stat});
          $currsheet->write($y, 4, $fixes->[1], $XLCOLORS{$stat});
          if ($BESTS{$pid} && $BESTS{$pid} eq $jid)
          {
            $bestsheet->write($y, 2, $loci, $XLCOLORS{$stat});
            $bestsheet->write($y, 3, $fixes->[0], $XLCOLORS{$stat});
            $bestsheet->write($y, 4, $fixes->[1], $XLCOLORS{$stat});
          }
        }
        $y++;
      }
      $x++;
    }
    $track->close();
    system "chmod 0777 $xlpath";
  }


}





print "brought to you by $JGISBV\n\n";

exit;

################################################################################
################################# SUBROUTINES ##################################
################################################################################

sub slurp
{
  my ($path) = @_;
  open (my $FH, '<', $path) || croak "Can't open $path: $OS_ERROR\n";
  my $raw = do {local $INPUT_RECORD_SEPARATOR = <$FH>};
  close $FH;
  my @arr = split m{[\n\r]}x, $raw;
  return @arr;
}

sub compare_variants
{
  my ($gvars, $ovars) = @_;

  my @gs = @{$gvars};
  my @os = @{$ovars};
  return () if (scalar @gs != scalar @os);
  my $lastpos = undef;
  my $farpos = undef;
  my $have = q{};
  while (scalar @os)
  {
    my $oline = shift @os;
    my $ohsh = variant_info($oline);

    my $gline = shift @gs;
    my $ghsh = variant_info($gline);

    carp ("\tbig oops, clones do not match") if ($ghsh->{clone} ne $ohsh->{clone});

    return () if ($ghsh->{clone} ne $ohsh->{clone});
    return () if ($ghsh->{offset} ne $ohsh->{offset});
    $lastpos = $lastpos || $ghsh->{offset};
    $farpos = $ghsh->{offset};
    return () if (abs($lastpos - $farpos) > 5);
    $have .= '(GATK: ' . $ghsh->{offset} . $ghsh->{want} . q{/};
    $have .= $ghsh->{have} . q{)};
    $have .= '(GCON: ' . $ohsh->{offset} . $ohsh->{want} . q{/};
    $have .= $ohsh->{have} . q{)};
  }

  return ($lastpos, $farpos, $have);
}

sub variant_info
{
  my ($line) = @_;
  my @atts = split m{\s}x, $line;
  my %hsh = (
    clone  => $atts[0],
    offset => $atts[1],
    want   => $atts[3],
    have   => $atts[4]
  );
  return \%hsh;
}

sub fixing_primers
{
  my ($fixarr, $refseq) = @_;
  my ($fstart, $fstop) = @{$fixarr};
  my $reflen = length $refseq;
  if (abs($reflen - $fstart) < $reflen/10)
  {
    return (undef, undef, undef);
  }
  my $bit = 1;
  my $uplen = 10;
  my $upstart = ($fstart - 1) - $uplen;
  my $upol = substr $refseq, $upstart, $uplen;
  my $upmelt = $GD->melt($upol);
  while ($upmelt < 58)
  {
    $upstart = $upstart - 1;
    $uplen++;
    $upol = substr $refseq, $upstart, $uplen;
    $upmelt = $GD->melt($upol);
  }
  while ($upmelt > 60.5)
  {
    $upstart++;
    $uplen = $uplen - 1;
    $upol = substr $refseq, $upstart, $uplen;
    $upmelt = $GD->melt($upol);
  }

  my $upcount = $GD->count($upol);
  $bit++ if $upcount->{GCp} > 58;
  $upol .= lc substr $refseq, $upstart + $uplen, 10;
  $upol = rcomp($upol);
  my $dnstart = $fstop;
  my $dnlen = 10;
  my $dnol = substr $refseq, $dnstart, $dnlen;
  my $dnmelt = $GD->melt($dnol);
  while ($dnmelt < 58)
  {
    $dnlen++;
    $dnol = substr $refseq, $dnstart, $dnlen;
    $dnmelt = $GD->melt($dnol);
  }
  while ($dnmelt > 60.5)
  {
    $dnlen = $dnlen - 1;
    $dnol = substr $refseq, $dnstart, $dnlen;
    $dnmelt = $GD->melt($dnol);
  }

  my $dncount = $GD->count($dnol);
  $bit++ if $dncount->{GCp} > 58;
  my $morednol = lc substr $refseq, $dnstart - 10, 10;
  $dnol = $morednol . $dnol;

  return ($bit, $upol, $dnol);
}

sub rcomp
{
  my ($strand) = @_;
  $strand = scalar reverse($strand);
  $strand =~ tr/AaCcGgTtRrYyKkMmBbDdHhVv/TtGgCcAaYyRrMmKkVvHhDdBb/;
  return $strand;
}

sub xmly
{
  my ($jid, $pid) = @_;
  my $prefix = substr $jid, 0, 3;
  my $gjid = $jid . '-GATK';
  my $locus = $pid;
  my $string = <<"END";
<?xml version="1.0" encoding="UTF-8" standalone="no"?>
<Session genome="http://synbio.jgi-psf.org/pacbio/analyses/splits/$pid.fasta" locus="$locus" version="4">
    <Resources>
      <Resource path="http://synbio.jgi-psf.org/pacbio/jobs/$prefix/$gjid/snps.ccs.vcf"/>
      <Resource path="http://synbio.jgi-psf.org/pacbio/jobs/$prefix/$gjid/realigned_reads.bam"/>
      <Resource path="http://synbio.jgi-psf.org/pacbio/jobs/$prefix/$jid/data/variants.bed"/>
      <Resource path="http://synbio.jgi-psf.org/pacbio/jobs/$prefix/$gjid/callable.bed"/>
      <Resource path="http://synbio.jgi-psf.org/pacbio/jobs/$prefix/$jid/data/aligned_reads.noref.bam"/>
    </Resources>
    <PanelLayout dividerFractions="0.06323529411764706,0.4235294117647059,0.8132352941176471"/>
</Session>

END
  return $string;
}



__END__

=head1 NAME

  JGISB_Analyze_PacBio.pl

=head1 VERSION

  Version 2.00

=head1 DESCRIPTION

  Decides which clones should be nominated from a PacBio run

=head1 USAGE

=head1 ARGUMENTS

Required arguments:

  -a,   --analysisids : A string of analysis numbers separated by commas.
              For example, -a 023175,023176,023177

Optional arguments:

  -h,   --help : Display this message

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2014, DNAssemble developers
All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright notice, this
list of conditions and the following disclaimer in the documentation and/or
other materials provided with the distribution.

* The names of the Joint Genome Institute, the Lawrence Berkeley National
Laboratory, the Department of Energy, and the JGI developers may not be used to
endorse or promote products derived from this software without specific prior
written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE DEVELOPERS BE LIABLE FOR ANY DIRECT, INDIRECT,
INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

=cut