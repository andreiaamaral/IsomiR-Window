#!/usr/bin/perl
#use warnings;
#use strict;
#use 5.010;
#
##This script aims to create an arf file to be used in the miRDeep2 meta-analysis project
##the first input is a fasta file, which resulted from the merging of all the collapsed fasta files of all the runs
## this file was created using the following command: home/local/opt/miRDeep/mirdeep2-0.1.0/bin/collapse_reads_md.pl <( cat *collapsed ) seq > all_colapsed_19Feb2019.fa
##then using bash I have created another file in TABLE format each line has the seq_id and the sequence
## the second input is file resulting from merging all the arf files, extracting the first column, sorting and selecting the unique lines
##create a merged arf file
##cat *_arf_file >> all_arf_merged &
##from this one remove first column, sort and select uniques
# #cat all_arf_merged | cut -f2,3,4,5,6,7,8,9,10,11,12,13,14 | sort | uniq > all_arf_merged_uniq & 
#  
#   #this script will create as output a all_arf_merged_uniq with the matching sequence ids from the all_colapsed_19Feb2019.fa
#    #created by AJAmaral Feb 2019 while at RTH
#     
my($path_to_all_collapsed_file, $path_to_all_arf_merged_uniq, $path_merged_arf_to_use_in_miRDeep2) = @ARGV;
my $collapsed_table = $ARGV[0];
my $merged_arf=$ARGV[1];
my $out_arf=$ARGV[2];
#open file with 2 columns: column1 seq id (index x number of reads) column2 sequence
open (GET_FILE_FASTA,"<$collapsed_table") or die "Can't open $collapsed_table: $!\n";
my %fasta;
#save content in hash key=sequence value = id
while (defined (my $line_table=<GET_FILE_FASTA>)){
	chomp $line_table;
	my @fields=split (' ', $line_table);
	my $id=$fields[0];
	my $sequence=$fields[1];
	my $sequence_lc= lc $sequence; 
	$fasta{$sequence_lc} = $id;
}
#open file with merged arf files, with only unique lines
open (GET_FILE_ARF, "<$merged_arf") or die "Can't open $merged_arf: $!\n";
open (OUT_ARF, ">$out_arf") or die "Can't write to $out_arf: $!\n";
#read through file and 
while (defined (my $line_arf=<GET_FILE_ARF>)){
	chomp $line_arf;
	my @fields2=split ('\t', $line_arf);
	my $seq=$fields2[3];
	if (exists $fasta{$seq}) {
	print "exists $seq\n";
	print OUT_ARF "$fasta{$seq}\t$line_arf\n";
	#print "$fasta{$seq}\t$line_arf\n";
	}
 }                                                                                                                                                   

close GET_FILE_FASTA;               	 									 	 	 	 	    
close GET_FILE_ARF;               	 									 	 	 	
close OUT_ARF;
exit;
