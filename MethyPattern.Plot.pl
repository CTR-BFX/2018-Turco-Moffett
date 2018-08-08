#!/usr/bin/perl
#------------------------------------------------------------------------------
# EPIC Analysis        
# Data from GenomeScan EPIC Array (organoids)
#
# Link to publication
# TO ADD ONCE AVAILABLE
#
# Script available from:
# https://github.com/CTR-BFX/2018-Turco-Moffett
#
#
# Analysis Performed by Russell S. Hamilton
# CTR Bioinformatics Facility
# Centre for Trophoblast Reseach, University of Cambridge, UK
# Copyright Russell S. Hamilton (rsh46@cam.ac.uk)
#
#------------------------------------------------------------------------------
# License Information
# This program is free software: you can redistribute it and/or modify  
# it under the terms of the GNU General Public License as published by  
# the Free Software Foundation, version 3.
#
# This program is distributed in the hope that it will be useful, but 
# WITHOUT ANY WARRANTY; without even the implied warranty of 
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License 
# along with this program. If not, see <http://www.gnu.org/licenses/>.
#------------------------------------------------------------------------------

use strict;

my $Gene         = "ELF5";
my $ROI          = "11:34500340-34535347"; # GRCh37 ELF5
my $Strand       = -1;
my $RefGenome    = "/usr/local/Genomes/Homo_sapiens/GRCh37.93/" . 
                   "Homo_sapiens.GRCh37.dna.chromosome.ALL.fa";
my $BedGraphLoc  = "Methylation_Files/";

my $Prom_Down = 0;
my $Prom_Up   = 400;

print $ROI, "\n";

my($ROI_CHR, $ROI_STA, $ROI_STO) = split(/[:-]/, $ROI);

my $strand_specific = "";
if($Strand == -1)
  { 
    print "Reverse Stranded\n";
    $ROI_STA = $ROI_STO - $Prom_Down;
    $ROI_STO = $ROI_STO + $Prom_Up;
  }

$ROI = "${ROI_CHR}:${ROI_STA}-${ROI_STO}";

print $ROI, "\n";

open(BED, ">tmp.bed") || die "Can't open tmp.bed for writing: $!\n";
print BED "$ROI_CHR\t$ROI_STA\t$ROI_STO\n";
close BED;

print "Extract Ref Seq\n";

my $ROI_Seq = `samtools faidx ${RefGenome} ${ROI} $strand_specific`;
chomp($ROI_Seq);
$ROI_Seq =~ s/^>.*//g;
$ROI_Seq =~ s/\n//g;



my $regex='CG';
my @CGLocations;
my %CGLocationsHash;
my @SampleCGs;

my $cgid = 1;
while($ROI_Seq =~ /($regex)/g) 
     {
       push @CGLocations, (pos($ROI_Seq)-length $1);
       $CGLocationsHash{(pos($ROI_Seq)-length $1)+$ROI_STA} = $cgid;
       $cgid++;
     }


print "+", "-"x79, "\n";
my $i; 
for($i=0; $i<=$#CGLocations; $i++)
   {
     printf "%3d  %5d  %10d %3d\n", $i+1, $CGLocations[$i], $CGLocations[$i]+$ROI_STA, 
                                    $CGLocationsHash{$CGLocations[$i]+$ROI_STA};
   }

print "+", "-"x79, "\n";

my @BedGraphFiles = ("org_1_BS.srt.bedgraph",  "org_3_BS.srt.bedgraph",
                     "org_4_BS.srt.bedgraph",  "org_5_BS.srt.bedgraph",
                     "pre_63_BS.srt.bedgraph", "pre_64_BS.srt.bedgraph",   
                     "pre_65_BS.srt.bedgraph", "pre_66_BS.srt.bedgraph");



for($i=0; $i<=$#BedGraphFiles; $i++)
   {
     print $BedGraphFiles[$i], "\n";
     my @intersect = `bedtools intersect -a ${BedGraphLoc}$BedGraphFiles[$i] -b tmp.bed`;

     foreach my $olap (@intersect)
       {
         chomp($olap);
         my($ol_chr, $ol_sta, $ol_sto, $ol_sco) = split(/\s+/, $olap);
         $SampleCGs[$i][$CGLocationsHash{$ol_sta}] = $ol_sco;
         printf "%3d  %10d  %6.4f\n", $CGLocationsHash{$ol_sta}, $ol_sta, $ol_sco;
       }
   }

print "+", "-"x79, "\n";

my $j;

print "CG";
for($i=0; $i<=$#BedGraphFiles; $i++){ print ", $BedGraphFiles[$i]"; }
print "\n";

for($i=0; $i<=$#CGLocations; $i++)
   {
     print $i+1;
  
     for($j=0; $j<=$#BedGraphFiles; $j++)
        {
          printf ", %6.4f", $SampleCGs[$j][$i];
        }
     print "\n";
   }


print "+", "-"x79, "\n";
print "+ END OF SCRIPT\n";
print "+", "-"x79, "\n";
