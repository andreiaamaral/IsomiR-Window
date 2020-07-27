#!/usr/bin/perl
use warnings;
use strict;
use Data::Dumper;
#
#
#################   DeSeq   ################################# 
#
## This script will receive the path to the folder that contains the .txt files 
## that resulted from running the script find_variation5.pl to every Filtered_SAM_file
## It will also receive the path to put the results and the information
## wether the samples are paired or unpaired
#
## These .txt files contain the miRNAs/isomiRs that were identified in each
## Filtered_SAM_file ant their respective frequency in the sample
#
## The script will run only one time
#
#
my($ID,$path_to_Deseq_C1_isomiRs, $path_to_Deseq_C2_isomiRs, $path_to_Deseq_C1_ncRNAs, $path_to_Deseq_C2_ncRNAs, $path_to_results, $paired_unpaired) = @ARGV; #receber mais dois para os ncRNAs
#
#
opendir my $dir_C1, "$path_to_Deseq_C1_isomiRs" or die "Cannot open directory: $!";
my @files_C1 = readdir $dir_C1;
closedir $dir_C1;

opendir my $dir_C2, "$path_to_Deseq_C2_isomiRs" or die "Cannot open directory: $!";
my @files_C2 = readdir $dir_C2;
closedir $dir_C2;
#
my $number_of_files_C1 = (scalar @files_C1) - 2;
my $number_of_files_C2 = (scalar @files_C2) - 2;
#print "@files_C1\n";
#print "@files_C2\n";
#splice @files_C1, $number_of_files_C1, 2;
#splice @files_C2, $number_of_files_C2, 2;
#foreach my $fileC1 (@files_C1){
#        chomp $fileC1;
#	print "$fileC1\n"
#}
my $i1 = 1;
my $i2 = 1;

print "after index OK\n";
 
my @names_files_C1 = ();
my @names_files_C2 = ();
my @header = ("IsomiR_ID");
my @design_C1 = ();
my @design_C2 = ();

foreach my $fileC1 (@files_C1){
	chomp $fileC1;
	if($fileC1 eq "." || $fileC1 eq ".."){ 
	next;
	}
	my $name_txt_file = $path_to_Deseq_C1_isomiRs."/".$fileC1;

	push(@names_files_C1, $name_txt_file);
	push(@header, "\tC1R$i1");
	push(@design_C1, "\nC1R$i1");
	$i1++;
}
print "1st foreach OK\n";
foreach my $fileC2 (@files_C2){
	chomp $fileC2;
        if($fileC2 eq "." || $fileC2 eq ".."){ 
        next;
        }
	my $name_txt_file = $path_to_Deseq_C2_isomiRs."/".$fileC2;	
	print "$name_txt_file\n";
	push(@names_files_C2, $name_txt_file);
	push(@header, "\tC2R$i2");
	push(@design_C2, "\nC2R$i2");
	$i2++;
}

print "2nd roreach OK\n";

push(@design_C1, @design_C2);

my @all_files = (@names_files_C1, @names_files_C2);

#print "@all_files\n";
my %files_hash;
my @array_frequecies = (0)x(scalar(@all_files));
my $num_of_file = 0;

foreach my $file (@all_files){

	open (GET_DATA_COUNTS,"<$file") or die "Can't open $file: $!\n";

	while (defined (my $line = <GET_DATA_COUNTS>)) {
		chomp $line;

		my ($isomiR_ID, $freq) = split("\t", $line);

		if (exists $files_hash{$isomiR_ID}){

			$files_hash{$isomiR_ID}[$num_of_file] = $freq;

		}

		else {

			$array_frequecies[$num_of_file] = $freq;
			$files_hash{$isomiR_ID} = [@array_frequecies];
		}

	}

	$num_of_file++;
}
#print Dumper(\%files_hash);
#
## Create the file with the counts of isomiRs for DESeq2 input
#
my $path_to_file = "$path_to_results"."/Deseq_IsomiR_$ID.txt";

open(OUTDESEQ,">$path_to_file") or die "Can't open $path_to_file: $!\n";

print OUTDESEQ "@header\n"; 

foreach my $key (keys %files_hash){

	my $string_of_array = join ( "\t", @{ $files_hash{$key} } );
	print OUTDESEQ "$key\t$string_of_array\n"; 

}

# Create the file with the experimental design for DESeq2 input for IsomiR analysis

my $path_experimental_design = "$path_to_results"."/ExperimentalDesign_$ID.IsomiR.txt";

if ($paired_unpaired eq "Unpaired_Sample") {

	open(OUTDESIGN,">$path_experimental_design") or die "Can't open $path_experimental_design: $!\n";

	foreach my $condition_replicate (@design_C1){

		my ($condition) = substr($condition_replicate, 0, index($condition_replicate, 'R'));
		$condition =~ s/\R//g;
		foreach my $filenameIso(@all_files){
			print OUTDESIGN "$condition_replicate\t$filenameIso\t$condition\n";
		}
	}

}

elsif ($paired_unpaired eq "Paired_Sample") {

	open(OUTDESIGN,">$path_experimental_design") or die "Can't open $path_experimental_design: $!\n";

	foreach my $condition_replicate (@design_C1){

		my ($condition) = substr($condition_replicate, 0, index($condition_replicate, 'R'));
		$condition =~ s/\R//g;
		my $subject = substr($condition_replicate, 3);
		foreach my $filenameIso(@all_files){																
		print OUTDESIGN "$condition_replicate\t$filenameIso\t$condition\t$subject\n";
		}
	}
															

}

#######################################################################################
#
## Do the exact same thing, but for the files of the other ncRNAs
#
opendir my $dir_C1_ncRNAs, "$path_to_Deseq_C1_ncRNAs" or die "Cannot open directory: $!";
my @files_C1_ncRNAs = readdir $dir_C1_ncRNAs;
closedir $dir_C1_ncRNAs;

opendir my $dir_C2_ncRNAs, "$path_to_Deseq_C2_ncRNAs" or die "Cannot open directory: $!";
my @files_C2_ncRNAs = readdir $dir_C2_ncRNAs;
closedir $dir_C2_ncRNAs;

my $number_of_files_C1_ncRNAs = (scalar @files_C1_ncRNAs) - 2;
my $number_of_files_C2_ncRNAs = (scalar @files_C2_ncRNAs) - 2;

#my $offsetC1=$number_of_files_C1_ncRNAs+1;
#splice @files_C1_ncRNAs, $number_of_files_C1_ncRNAs, 2;
#splice @files_C2_ncRNAs, $number_of_files_C2_ncRNAs, 2;


#
my $i3 = 1;
my $i4 = 1;
#
#
my @names_files_C1_ncRNAs = ();
my @names_files_C2_ncRNAs = ();
my @header_ncRNAs = ("ncRNA_ID");
#
print "after index OK\n";


foreach my $fileC1_ncRNA (@files_C1_ncRNAs){
        chomp $fileC1_ncRNA;
	if($fileC1_ncRNA eq "." || $fileC1_ncRNA eq ".."){
        next;
        }
	my $name_txt_file = $path_to_Deseq_C1_ncRNAs."/".$fileC1_ncRNA;
	push(@names_files_C1_ncRNAs, $name_txt_file);
        push(@header_ncRNAs, "\tC1R$i3");
        push(@design_C1, "\nC1R$i3");
        $i3++;
}
print "ncRNA 1st foreach OK\n";

foreach my $fileC2_ncRNA (@files_C2_ncRNAs){
	chomp $fileC2_ncRNA;
	if($fileC2_ncRNA eq "." || $fileC2_ncRNA eq ".."){
        next;
        }
        my $name_txt_file = $path_to_Deseq_C2_ncRNAs."/".$fileC2_ncRNA;
        push(@names_files_C2_ncRNAs, $name_txt_file);
        push(@header_ncRNAs, "\tC2R$i4");
        push(@design_C2, "\nC2R$i4");
        $i4++;
}																			
#
my @all_files_ncRNAs = (@names_files_C1_ncRNAs, @names_files_C2_ncRNAs);

my %files_hash_ncRNAs;
my @array_frequecies_ncRNAs = (0)x(scalar(@all_files_ncRNAs));
my $num_of_file_ncRNAs = 0;

foreach my $file_ncRNAs (@all_files_ncRNAs){

	open (GET_DATA_COUNTS_NCRNAS,"<$file_ncRNAs") or die "Can't open $file_ncRNAs: $!\n";

	while (defined (my $line = <GET_DATA_COUNTS_NCRNAS>)) {
		chomp $line;
		my ($ncRNA_ID, $freq) = split("\t", $line);																				  if (exists $files_hash_ncRNAs{$ncRNA_ID}){
                    $files_hash_ncRNAs{$ncRNA_ID}[$num_of_file_ncRNAs] = $freq;
     																										          }
	        else {
																											        $array_frequecies_ncRNAs[$num_of_file_ncRNAs] = $freq;																			  $files_hash_ncRNAs{$ncRNA_ID} = [@array_frequecies_ncRNAs];
																										          }
         }
         $num_of_file_ncRNAs++;
}																						
#print Dumper(\%files__hash);
# Create the file with the counts of ncRNAs for DESeq2 input 																	

my $path_to_file_ncRNAs = "$path_to_results"."/Deseq_ncRNAs_$ID.txt";
																	 open(OUTDESEQ_NC,">$path_to_file_ncRNAs") or die "Can't open $path_to_file_ncRNAs: $!\n";
																									
print OUTDESEQ_NC "@header_ncRNAs\n";

																	 foreach my $key_ncRNAs (keys %files_hash_ncRNAs){
	my $string_of_array_ncRNAs = join ( "\t", @{ $files_hash_ncRNAs{$key_ncRNAs} } );
	print OUTDESEQ_NC  "$key_ncRNAs\t$string_of_array_ncRNAs\n"; 

}

# Create the file with the experimental design for DESeq2 input for ncRNA analysis
#
 my $path_experimental_design2 = "$path_to_results"."/ExperimentalDesign_$ID.ncRNAs.txt";
#
 if ($paired_unpaired eq "Unpaired_Sample") {
#
         open(OUTDESIGN2,">$path_experimental_design2") or die "Can't open $path_experimental_design: $!\n";
#
	foreach my $condition_replicate (@design_C1){

              my ($condition) = substr($condition_replicate, 0, index($condition_replicate, 'R'));
              $condition =~ s/\R//g;
              foreach my $filenameNC(@all_files_ncRNAs){
                     print OUTDESIGN2 "$condition_replicate\t$filenameNC\t$condition\n";
               }
         }

}

elsif ($paired_unpaired eq "Paired_Sample") {

        open(OUTDESIGN2,">$path_experimental_design2") or die "Can't open $path_experimental_design: $!\n";

        foreach my $condition_replicate (@design_C1){

                my ($condition) = substr($condition_replicate, 0, index($condition_replicate, 'R'));
                $condition =~ s/\R//g;
                my $subject = substr($condition_replicate, 3);
                foreach my $filenameNC(@all_files_ncRNAs){                                                                                     
                print OUTDESIGN2 "$condition_replicate\t$filenameNC\t$condition\t$subject\n";
                }
        }


}

#
exit;
