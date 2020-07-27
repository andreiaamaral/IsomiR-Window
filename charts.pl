#!/usr/bin/perl
use warnings;
use strict;
use Data::Dumper;
use List::Util qw[min max];
use POSIX;

require "./Species.pl";

my($ID, $path_to_Deseq_C1_isomiRs, $path_to_Deseq_C2_isomiRs, $path_to_Deseq_C1_isomiRs_amb, $path_to_Deseq_C2_isomiRs_amb, 
	$path_to_ncRNAs_C1, $path_to_ncRNAs_C2, $path_to_files_charts_C1, $path_to_files_charts_C2, $path_to_general_results, $species) = @ARGV;


#########################################################################################################################
# Get the number of not_ambiguous.txt files in the specified folder
my $number_of_files_C1 = get_number_files_directory($path_to_Deseq_C1_isomiRs);

my $number_of_files_C2 = get_number_files_directory($path_to_Deseq_C2_isomiRs);

# Inititate counter variables
my $i1 = 1;
my $i2 = 1;

# Create arrays to store the names of the files for each condition
my @names_files_C1 = ();
my @names_files_C2 = ();


# Go trhough the folders and store the names of the files in an array for each condition
while ($i1 <= $number_of_files_C1){

	my $name_txt_file = "$path_to_Deseq_C1_isomiRs/Filtered_SAM_"."$ID"."_C1_R"."$i1".".fastq.sam_not_ambiguous.txt";

	push(@names_files_C1, $name_txt_file);
	$i1++;
}

while ($i2 <= $number_of_files_C2){

	my $name_txt_file = "$path_to_Deseq_C2_isomiRs/Filtered_SAM_"."$ID"."_C2_R"."$i2".".fastq.sam_not_ambiguous.txt";

	push(@names_files_C2, $name_txt_file);
	$i2++;
}

my %hash_total_isomiRs_C1; 
my %hash_total_isomiRs_C2; 
print %hash_total_isomiRs_C1;
print %hash_total_isomiRs_C2;

%hash_total_isomiRs_C1 = create_hash_counts("pie", @names_files_C1, %hash_total_isomiRs_C1);
%hash_total_isomiRs_C2 = create_hash_counts("pie", @names_files_C2, %hash_total_isomiRs_C2);

#print Dumper(\%hash_total_isomiRs_C1);
#print Dumper(\%hash_total_isomiRs_C2);

my $total_num_miRNAs_C1 = 0;
my $total_num_miRNAs_C2 = 0;

$total_num_miRNAs_C1 = add_miRNAs($total_num_miRNAs_C1, %hash_total_isomiRs_C1);
$total_num_miRNAs_C2 = add_miRNAs($total_num_miRNAs_C2, %hash_total_isomiRs_C2);


#print "$total_num_miRNAs_C1\n";
#print "$total_num_miRNAs_C2\n";

#print Dumper(\%hash_total_isomiRs_C1);


##### Do the same procedure to the ambiguous miRNAs files and add the counts to the total miRNAs variable

# Get the number of ambiguous.txt files in the specified folder
my $number_of_files_C1_isomiRs_amb = get_number_files_directory($path_to_Deseq_C1_isomiRs_amb);
my $number_of_files_C2_isomiRs_amb = get_number_files_directory($path_to_Deseq_C2_isomiRs_amb);

# Inititate counter variables
my $i3 = 1;
my $i4 = 1;

# Create arrays to store the names of the files for each condition
my @names_files_C1_isomiRs_amb = ();
my @names_files_C2_isomiRs_amb = ();


# Go through the folders and store the names of the files in an array for each condition
while ($i3 <= $number_of_files_C1_isomiRs_amb){

	my $name_txt_file = "$path_to_Deseq_C1_isomiRs_amb/Filtered_SAM_"."$ID"."_C1_R"."$i3".".fastq.sam_ambiguous.txt";

	push(@names_files_C1_isomiRs_amb, $name_txt_file);
	$i3++;
}

while ($i4 <= $number_of_files_C2_isomiRs_amb){

	my $name_txt_file = "$path_to_Deseq_C2_isomiRs_amb/Filtered_SAM_"."$ID"."_C2_R"."$i4".".fastq.sam_ambiguous.txt";

	push(@names_files_C2_isomiRs_amb, $name_txt_file);
	$i4++;
}

# Define the hashes to store all the miRNAs throughout the files and sum the frequencies
my %hash_total_isomiRs_C1_isomiRs_amb;
my %hash_total_isomiRs_C2_isomiRs_amb;

%hash_total_isomiRs_C1_isomiRs_amb = create_hash_counts("pie", @names_files_C1_isomiRs_amb, %hash_total_isomiRs_C1_isomiRs_amb);
%hash_total_isomiRs_C2_isomiRs_amb = create_hash_counts("pie", @names_files_C2_isomiRs_amb, %hash_total_isomiRs_C2_isomiRs_amb);

#print Dumper(\%hash_total_isomiRs_C1_isomiRs_amb);
#print Dumper(\%hash_total_isomiRs_C2_isomiRs_amb);

# Go through the created hashes and add the ambiguous miRNAs to the total 
# number of miRNAs for each condition 
$total_num_miRNAs_C1 = add_miRNAs($total_num_miRNAs_C1, %hash_total_isomiRs_C1_isomiRs_amb);
$total_num_miRNAs_C2 = add_miRNAs($total_num_miRNAs_C2, %hash_total_isomiRs_C2_isomiRs_amb);

########################## Produce the files to create the pie chart of the total number of ncRNAs ######################

# Get the number of ncRNAs.txt files in the specified folder 

my $number_of_files_C1_ncRNAs = get_number_files_directory($path_to_ncRNAs_C1);
my $number_of_files_C2_ncRNAs = get_number_files_directory($path_to_ncRNAs_C2);

# Inititate counter variables
my $i5 = 1;
my $i6 = 1;

# Create arrays to store the names of the files for each condition
my @names_files_C1_ncRNAs = ();
my @names_files_C2_ncRNAs = ();

# Go trhough the folders and store the names of the files in an array for each condition
while ($i5 <= $number_of_files_C1_ncRNAs){

	my $name_txt_file = "$path_to_ncRNAs_C1/Counts_ncRNAs_HTSeq_"."$ID"."_C1_R"."$i5".".fastq.txt";

	push(@names_files_C1_ncRNAs, $name_txt_file);
	$i5++;
}

while ($i6 <= $number_of_files_C2_ncRNAs){

	my $name_txt_file = "$path_to_ncRNAs_C2/Counts_ncRNAs_HTSeq_"."$ID"."_C2_R"."$i6".".fastq.txt";

	push(@names_files_C2_ncRNAs, $name_txt_file);
	$i6++;
}

# Define the hashes to store all the ncRNAs trhoughout the files and sum the frequencies
my %hash_total_isomiRs_C1_ncRNAs;
my %hash_total_isomiRs_C2_ncRNAs;

%hash_total_isomiRs_C1_ncRNAs = create_hash_counts("pie", @names_files_C1_ncRNAs, %hash_total_isomiRs_C1_ncRNAs);
%hash_total_isomiRs_C2_ncRNAs = create_hash_counts("pie", @names_files_C2_ncRNAs, %hash_total_isomiRs_C2_ncRNAs);


#print Dumper(\%hash_total_isomiRs_C1_ncRNAs);
#print Dumper(\%hash_total_isomiRs_C2_ncRNAs);

my $path_output_total_ncRNAs_C1 = "$path_to_files_charts_C1/Total_ncRNAs_piechart_C1.txt";
my $path_output_total_ncRNAs_C2 = "$path_to_files_charts_C2/Total_ncRNAs_piechart_C2.txt";

# Create the output file to construct the pie chart for each condition
open(OUTDESEQ_NC_C1,">$path_output_total_ncRNAs_C1") or die "Can't open $path_output_total_ncRNAs_C1: $!\n";
open(OUTDESEQ_NC_C2,">$path_output_total_ncRNAs_C2") or die "Can't open $path_output_total_ncRNAs_C2: $!\n";

# Create an array of all possible ncRNAs in the samples
my @all_categories_ncRNAs_C1 = ("antisense_RNA", "guide_RNA", "lncRNA", "misc_RNA", "other", "ribozyme", "RNase_MRP_RNA", "RNase_P_RNA", 
	"rRNA", "scRNA", "snoRNA", "snRNA", "SRP_RNA", "telomerase_RNA", "tRNA", "vault_RNA");

# Write the counts to the file of Condition 1
foreach my $key5 (sort keys %hash_total_isomiRs_C1_ncRNAs){

	if (grep( $_ eq $key5, @all_categories_ncRNAs_C1)) {

		print OUTDESEQ_NC_C1 "$key5".","."$hash_total_isomiRs_C1_ncRNAs{$key5}\n";

		@all_categories_ncRNAs_C1 = grep { $_ ne "$key5" } @all_categories_ncRNAs_C1;
	}
}

foreach (@all_categories_ncRNAs_C1){

	print OUTDESEQ_NC_C1 "$_".","."0\n";
}

print OUTDESEQ_NC_C1 "miRNA".","."$total_num_miRNAs_C1\n";

# Create an array of all possible ncRNAs in the samples
my @all_categories_ncRNAs_C2 = ("antisense_RNA", "guide_RNA", "lncRNA", "misc_RNA", "other", "ribozyme", "RNase_MRP_RNA", "RNase_P_RNA", 
	"rRNA", "scRNA", "snoRNA", "snRNA", "SRP_RNA", "telomerase_RNA", "tRNA", "vault_RNA");

# Write the counts to the file of Condition 2
foreach my $key6 (sort keys %hash_total_isomiRs_C2_ncRNAs){

	if (grep( $_ eq $key6, @all_categories_ncRNAs_C2)) {

		print OUTDESEQ_NC_C2 "$key6".","."$hash_total_isomiRs_C2_ncRNAs{$key6}\n";

		@all_categories_ncRNAs_C2 = grep { $_ ne "$key6" } @all_categories_ncRNAs_C2;
	}
}


foreach (@all_categories_ncRNAs_C2){

	print OUTDESEQ_NC_C2 "$_".","."0\n";
}


print OUTDESEQ_NC_C2 "miRNA".","."$total_num_miRNAs_C2\n";

########################## Produce the files to create the pie chart of isomiRs variation ######################


my $path_output_types_isomirs_C1 = "$path_to_files_charts_C1/Types_isomiRs_C1.txt";
my $path_output_types_isomirs_C2 = "$path_to_files_charts_C2/Types_isomiRs_C2.txt";

my @types_isomirs = ("5' end addition", "5' end trimming", "3' end addition", "3' end trimming", "5' end addition and 3' end trimming", "5' end trimming and 3' end addition", "5' end and 3' end trimming", "5' end and 3' end addition", "internal editings", "tailings");

my @counts_types_isomirs_C1 = count_types_isomirs(%hash_total_isomiRs_C1);
my @counts_types_isomirs_C2 = count_types_isomirs(%hash_total_isomiRs_C2);

open(OUT_TYPES_ISOMIRS_C1, ">$path_output_types_isomirs_C1") or die "Can't open $path_output_types_isomirs_C1: $!\n";
open(OUT_TYPES_ISOMIRS_C2, ">$path_output_types_isomirs_C2") or die "Can't open $path_output_types_isomirs_C2: $!\n";

my $counter = 0;

foreach my $type (@types_isomirs){

	print OUT_TYPES_ISOMIRS_C1 "$type".","."$counts_types_isomirs_C1[$counter]\n";
	print OUT_TYPES_ISOMIRS_C2 "$type".","."$counts_types_isomirs_C2[$counter]\n";

	$counter++;
}

########################## Produce the files to create the pie chart of 5p/3p ######################

my $path_output_3p5p_C1 = "$path_to_files_charts_C1/Counts_3p5p_C1.txt";
my $path_output_3p5p_C2 = "$path_to_files_charts_C2/Counts_3p5p_C2.txt";

my @counts_miRNAs_3p5p_C1 = count_3p5p(%hash_total_isomiRs_C1);
my @counts_miRNAs_3p5p_C2 = count_3p5p(%hash_total_isomiRs_C2);

open(OUT_ISOMIRS_3P5P_C1, ">$path_output_3p5p_C1") or die "Can't open $path_output_3p5p_C1: $!\n";
open(OUT_ISOMIRS_3P5P_C2, ">$path_output_3p5p_C2") or die "Can't open $path_output_3p5p_C2: $!\n";

my @p3p5_isomirs = ("3p arm", "5p arm");

my $counter_3p5p = 0;

foreach my $type_p3p5 (@p3p5_isomirs){

	print OUT_ISOMIRS_3P5P_C1 "$type_p3p5".","."$counts_miRNAs_3p5p_C1[$counter_3p5p]\n";
	print OUT_ISOMIRS_3P5P_C2 "$type_p3p5".","."$counts_miRNAs_3p5p_C2[$counter_3p5p]\n";

	$counter_3p5p++;
}

########################## Produce the files to create the bar chart of the types of editings ######################

my %hash_isomiRs_C1; 
my %hash_isomiRs_C2; 

%hash_isomiRs_C1 = create_hash_counts("bar", @names_files_C1, %hash_isomiRs_C1);
%hash_isomiRs_C2 = create_hash_counts("bar", @names_files_C2, %hash_isomiRs_C2);


my $path_output_editings_C1 = "$path_to_files_charts_C1/Counts_editings_C1.txt";
my $path_output_editings_C2 = "$path_to_files_charts_C2/Counts_editings_C2.txt";

my %counts_miRNAs_editings_raw_C1 = count_editings(%hash_isomiRs_C1);
my %counts_miRNAs_editings_raw_C2 = count_editings(%hash_isomiRs_C2);

# Calculate average and standard deviation
my %counts_miRNAs_editings_C1 = calculate_stDev(%counts_miRNAs_editings_raw_C1);
my %counts_miRNAs_editings_C2 = calculate_stDev(%counts_miRNAs_editings_raw_C2);

open(OUT_ISOMIRS_ED_C1, ">$path_output_editings_C1") or die "Can't open $path_output_editings_C1: $!\n";
open(OUT_ISOMIRS_ED_C2, ">$path_output_editings_C2") or die "Can't open $path_output_editings_C2: $!\n";


foreach my $type_editings (keys %counts_miRNAs_editings_C1){

	print OUT_ISOMIRS_ED_C1 "$type_editings".","."$counts_miRNAs_editings_C1{$type_editings}\n";
	print OUT_ISOMIRS_ED_C2 "$type_editings".","."$counts_miRNAs_editings_C2{$type_editings}\n";

}

########################## Produce the files to create the bar chart of the types of tailings ######################

my $path_output_tailings_C1 = "$path_to_files_charts_C1/Counts_tailings_C1.txt";
my $path_output_tailings_C2 = "$path_to_files_charts_C2/Counts_tailings_C2.txt";

my %counts_miRNAs_tailings_raw_C1 = count_tailings(%hash_isomiRs_C1);
my %counts_miRNAs_tailings_raw_C2 = count_tailings(%hash_isomiRs_C2);

# Calculate average and standard deviation
my %counts_miRNAs_tailings_C1 = calculate_stDev(%counts_miRNAs_tailings_raw_C1);
my %counts_miRNAs_tailings_C2 = calculate_stDev(%counts_miRNAs_tailings_raw_C2);

open(OUT_ISOMIRS_TAIL_C1, ">$path_output_tailings_C1") or die "Can't open $path_output_tailings_C1: $!\n";
open(OUT_ISOMIRS_TAIL_C2, ">$path_output_tailings_C2") or die "Can't open $path_output_tailings_C2: $!\n";


foreach my $type_tailings (keys %counts_miRNAs_tailings_C1){

	print OUT_ISOMIRS_TAIL_C1 "$type_tailings".","."$counts_miRNAs_tailings_C1{$type_tailings}\n";
	print OUT_ISOMIRS_TAIL_C2 "$type_tailings".","."$counts_miRNAs_tailings_C2{$type_tailings}\n";

}

########################## Produce the files to create the final DeSeq table and the bar chart of the miRNAs families ######################

my $path_output_families_C1 = "$path_to_files_charts_C1/Counts_families_C1.txt";
my $path_output_families_C2 = "$path_to_files_charts_C2/Counts_families_C2.txt";


my $table_de = $path_to_general_results."/TableDE_".$ID.".txt";
unless (-e $table_de) {
	$table_de = $path_to_general_results."/CompleteTableDE_".$ID.".txt";
}
my $final_table = $path_to_general_results."/Final_TableDE_".$ID.".txt";
my $miRNA_normalized_counts_file = $path_to_general_results."/miRNA_normalized_counts_".$ID.".txt";

open(TABLE_DE, '<', $table_de) or die "Could not open $table_de: $!";

# Get the header
my $de_header = <TABLE_DE>;

chomp $de_header;
my @header_parts = split(/\t/, $de_header);
#print @header_parts;
my $number_C1_files;
my $number_C2_files;
my $num_columns;

foreach my $header_column (@header_parts) {
	$num_columns = scalar @header_parts;

	my $replicate = substr $header_column,0,3;
	#print $replicate;
	if ($replicate eq "C1R") {
		$number_C1_files++;
	}
	elsif ($replicate eq "C2R") {
		$number_C2_files++;
	}
}

my @lines_with_average_stdev;

my %miRNAs_by_replicate_C1;
my %miRNAs_by_replicate_C2;

my %miRNAs_by_replicate;

# Create a new DE table with the average and stdev for each condition, for each miRNA
foreach my $de_line (<TABLE_DE>) {
		chomp $de_line;
		my @de_row = split(/\t/, $de_line);
		
		my @isomiR_id = split('_',$de_row[0]);
		my $miRNA = $isomiR_id[1];

		my @C1_values;
		my @C2_values;

		# Go through each column of the line and fill the hash with the miRNA name as key and
		# the sum of all miRNAs/isomiRs of that type for each replicate (column)
		# Ex: a file with 2 replicates for each condition, each key/value pair would look like:
		# hsa-miR-183-5p : [22, 41, 35, 36] (if C1_R1 has a sum of 22 counts, C1_R2 has a sum of 41 counts, etc)
		my @miRNA_values;
		my $dict_position;
		if (exists $miRNAs_by_replicate{$miRNA}) {
			$dict_position = 0;
			foreach my $replicate_column (7 .. $num_columns-1) {
				@{$miRNAs_by_replicate{$miRNA}}[$dict_position] += $de_row[$replicate_column];
				$dict_position++;
			}
		}
		else {
			$dict_position = 0;
			foreach my $replicate_column (7 .. $num_columns-1) {

				@{$miRNAs_by_replicate{$miRNA}}[$dict_position] = $de_row[$replicate_column];
				$dict_position++;
			}


		}

		# Get the counts by isomiR, by experimental condition
		foreach my $C1_replicate (0 .. $number_C1_files-1) {
			# The replicate information for C1 starts after the 7th line
			my $C1_value = $de_row[7+$C1_replicate];
			push(@C1_values, $C1_value);

		}


		foreach my $C2_replicate (0 .. $number_C2_files-1) {
			# The replicate information for C2 starts after the 7th line + number of C1 replicates
			my $C2_value = $de_row[7+$number_C1_files+$C2_replicate];
			push(@C2_values, $C2_value);

		}

		# Get the average and stdev for each line of each condition
		my $C1_average = average(@C1_values);
		my $C1_stdev = stdev(@C1_values);

		my $C2_average = average(@C2_values);
		my $C2_stdev = stdev(@C2_values);

		splice @de_row, 7;
		push(@de_row, $C1_average, $C1_stdev, $C2_average, $C2_stdev);
		my $final_line = join("\t", @de_row);
		push(@lines_with_average_stdev, $final_line);
}

open(FINAL_TABLE_DE, '>', $final_table) or die "Could not open $final_table: $!";
my $final_header_table = join("\t", "IsomiRs","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj","C1_average","C1_stDev","C2_average","C2_stDev");
print FINAL_TABLE_DE "$final_header_table\n";

my %hash_DE_C1;
my %hash_DE_C2;

# Write the average and stDev of each isomiR for each experimental condition in a new file
foreach my $final_header_line (@lines_with_average_stdev) {
	print FINAL_TABLE_DE "$final_header_line\n";
}

# Write the miRNAs normalized counts file
open(MIRNA_NORMALIZED_COUNTS, '>', $miRNA_normalized_counts_file) or die "Could not open $miRNA_normalized_counts_file: $!";

my @normalized_header = "miRNA_name";
foreach my $c1files(1 .. $number_C1_files) { push @normalized_header, "C1R".$c1files; }
foreach my $c2files (1 .. $number_C2_files) { push @normalized_header, "C2R".$c2files; }
my $final_norm_header = join("\t",@normalized_header);
print MIRNA_NORMALIZED_COUNTS $final_norm_header, "\n";

foreach my $miRNA_replicate_line (sort keys %miRNAs_by_replicate) {
	my @rep_values_list = @{$miRNAs_by_replicate{$miRNA_replicate_line}};
	my $rep_values_str = join("\t", @rep_values_list);

	print MIRNA_NORMALIZED_COUNTS $miRNA_replicate_line,"\t", $rep_values_str, "\n";
}

close(MIRNA_NORMALIZED_COUNTS);

# For each condition, for each miRNA family calculate the average counts and stdev and write it to a file
open(FAMILIES_C1, '>', $path_output_families_C1) or die "Could not open $path_output_families_C1: $!";
open(FAMILIES_C2, '>', $path_output_families_C2) or die "Could not open $path_output_families_C2: $!";

foreach my $family_mirna (sort keys %miRNAs_by_replicate) {
	#print($miRNAs_by_replicate{$family_mirna});
	my @C2_replicates = splice @{$miRNAs_by_replicate{$family_mirna}}, $number_C1_files;

	my @C1_replicates = @{$miRNAs_by_replicate{$family_mirna}}; #rest of the List
	my $C1_fam_avg = average(@C1_replicates);
	my $C1_fam_stdev = stdev(@C1_replicates);
	my $fam_C1_line = $family_mirna.','.$C1_fam_avg.'_'.$C1_fam_stdev."\n";
	print FAMILIES_C1 $fam_C1_line;

	my $C2_fam_avg = average(@C2_replicates);
	my $C2_fam_stdev = stdev(@C2_replicates);
	my $fam_C2_line = $family_mirna.','.$C2_fam_avg.'_'.$C2_fam_stdev."\n";
	print FAMILIES_C2 $fam_C2_line;
	
}
close(FAMILIES_C1);
close(FAMILIES_C2);

close(TABLE_DE);
close(FINAL_TABLE_DE);

########################## Produce the files to create the bar chart of the read size ##################################

# Inititate counter variables
my $i7 = 1;
my $i8 = 1;

# Create arrays to store the names of the files for each condition
my @names_files_C1_alignment = ();
my @names_files_C2_alignment = ();

# Go through the folders and store the names of the files in an array for each condition
while ($i7 <= $number_of_files_C1_ncRNAs){

	my $name_txt_file = "$path_to_files_charts_C1/"."$ID"."_C1_R"."$i7"."/"."$ID"."_C1_R"."$i7".".fastq.sam";

	push(@names_files_C1_alignment, $name_txt_file);
	$i7++;
}

while ($i8 <= $number_of_files_C2_ncRNAs){

	my $name_txt_file = "$path_to_files_charts_C2/"."$ID"."_C2_R"."$i8"."/"."$ID"."_C2_R"."$i8".".fastq.sam";

	push(@names_files_C2_alignment, $name_txt_file);
	$i8++;
}

# Define the hashes to store all the sizes of reads throughout the alignment files and sum the frequencies
my %hash_read_size_C1_raw = read_sizes(@names_files_C1_alignment);
my %hash_read_size_C2_raw = read_sizes(@names_files_C2_alignment);

my %hash_read_size_C1 = calculate_stDev(%hash_read_size_C1_raw);
my %hash_read_size_C2 = calculate_stDev(%hash_read_size_C2_raw);
#print Dumper(\%hash_read_size_C1);

#print Dumper(\%hash_read_size_C2);

my $path_output_read_size_C1 = "$path_to_files_charts_C1/Read_size_barchart_C1.txt";
my $path_output_read_size_C2 = "$path_to_files_charts_C2/Read_size_barchart_C2.txt";

# Create the output file to construct the bar chart for each condition
open(OUT_READSIZE_C1,">$path_output_read_size_C1") or die "Can't open $path_output_read_size_C1: $!\n";
open(OUT_READSIZE_C2,">$path_output_read_size_C2") or die "Can't open $path_output_read_size_C2: $!\n";


# Sorts the hash by the average values and writes the read length + ',' average_stdDev into a file
foreach my $readSize_C1 (sort keys %hash_read_size_C1) {
    print  OUT_READSIZE_C1 "$readSize_C1".","."$hash_read_size_C1{$readSize_C1}\n";
}

# Sorts the hash by the average values and writes the read length + ',' average_stdDev into a file
foreach my $readSize_C2 (sort keys %hash_read_size_C2) {
    print  OUT_READSIZE_C2 "$readSize_C2".","."$hash_read_size_C2{$readSize_C2}\n";
}

close OUT_READSIZE_C1;
close OUT_READSIZE_C2;

####################### Produce the files to create the isomiR positions chart ###########################
my $path_output_positions_C1 = "$path_to_files_charts_C1/isomiR_positions_C1.txt";
my $path_output_positions_C2 = "$path_to_files_charts_C2/isomiR_positions_C2.txt";

my $kingdom = Species::species_to_column($species, "kingdom");
my $positions1_maxVal = calculate_positions($path_output_positions_C1,$kingdom,@names_files_C1);
my $positions2_maxVal = calculate_positions($path_output_positions_C2,$kingdom,@names_files_C2);
my $positions_maxVal = max($positions1_maxVal,$positions2_maxVal);


########################## Produce the file to set the scale of the bar charts ###########################

my @readSize_files = ($path_output_read_size_C1, $path_output_read_size_C2);
my @families_files = ($path_output_families_C1, $path_output_families_C2);


my $readSize_maxVal = get_scales(@readSize_files);
my $families_maxVal = get_scales(@families_files);


my $file_chart_scales = $path_to_general_results . '/Chart_scales_' . $ID. '.txt';

# Create the output file with the max scales for each pair of charts
open(OUT_SCALES,">$file_chart_scales") or die "Can't open $file_chart_scales: $!\n";
print OUT_SCALES "Read_sizes:" . "$readSize_maxVal\n";
print OUT_SCALES "miRNA_families:" . "$families_maxVal\n";
print OUT_SCALES "isomiR_positions:" . "$positions_maxVal\n";

close OUT_SCALES;

############################################### SUBROUTINES ##############################################

sub get_number_files_directory {

	# This subroutine receives a path to a directory and returns the number of files
	# that the folder contains

	my($path) = (@_);

	opendir my $dir, "$path" or die "Cannot open directory: $!";
	my @files = readdir $dir;
	closedir $dir;

	my $number_of_files = (scalar @files) - 2;

	return $number_of_files;
}


sub create_hash_counts {

	# This subroutine goes through an array of files,
	# opens each file and creates a hash in which the first
	# column is the key and the second is the value

	# if the chart type is pie, returns a hash like 'element => average of values'
	# if the chart type is bar, return a hash like 'element => value_value_value'
	my ($chart_type, @array_names, %hash_counts) = (@_);
	my $num_of_replicates = scalar @array_names;


	foreach my $file (@array_names) {

		open (GET_DATA_COUNTS, "<$file") or die "Can't open $file: $!\n";

		while (defined (my $line = <GET_DATA_COUNTS>)) {

			chomp $line;

			my (@line_content) = split("\t", $line);

			my $element = $line_content[0];
			my $freq = $line_content[1];

			if (exists $hash_counts{$element}){

				$hash_counts{$element} = $hash_counts{$element} . '_' . $freq;
			}

			else {

				$hash_counts{$element} = $freq;
			}

		}
	}

	foreach my $key (keys %hash_counts){
		if ($chart_type eq 'pie') {
			my $value_sum;
			my @hash_values = split('_', $hash_counts{$key});

			foreach my $hash_value (@hash_values) {
				$value_sum += $hash_value;
			}

			# get the average of the values
			$hash_counts{$key} = $value_sum/$num_of_replicates;
		}
	}

	return %hash_counts;
}

sub add_miRNAs {

	# This subroutine receives a counter variable and a hash
	# Goes through the hash and adds each value of the key to 
	# the counter variable

	my ($counter, %hash_miRNAs) = (@_);

	foreach my $key (keys %hash_miRNAs) {

		$counter = $counter + $hash_miRNAs{$key};
	}

	return $counter;

}

sub count_types_isomirs {

	# This subroutine receives a hash of counts and returns an array with the
	# the total counts for each type of isomiRs

	my (%hash_isomirs) = (@_);

	my $a5 = 0;
	my $t5 = 0;
	my $a3 = 0;
	my $t3 = 0;
	my $a5t3 = 0;
	my $t5a3 = 0;
	my $t5t3 = 0;
	my $a5a3 = 0;
	my $editings = 0;
	my $tailings = 0;

	foreach my $key_isomir (keys %hash_isomirs) {


		my @info_id = split("_", $key_isomir);

		my $multiple_precursors = $info_id[2];

		if ($multiple_precursors ne 'multiple-precursors') {

			# If it has 3' additions
			if ($info_id[2] ne "0") {

				# Search for the isomiRs that only have 3' additions 
				if ($info_id[3] eq "0" and $info_id[4] eq "0" and $info_id[5] eq "0" and $info_id[6] eq "NA" and $info_id[9] eq "F"){

					$a3 = $a3 + $hash_isomirs{$key_isomir};

				}

				# Search for the isomiRs that have both 3' additions and 5' additions
				elsif ($info_id[3] eq "0" and $info_id[4] ne "0" and $info_id[5] eq "0" and $info_id[6] eq "NA" and $info_id[9] eq "F"){

					$a5a3 = $a5a3 + $hash_isomirs{$key_isomir};

				}

				# Search for the isomiRs that have both 3' additions and 5' trimmings
				elsif ($info_id[3] eq "0" and $info_id[4] eq "0" and $info_id[5] ne "0" and $info_id[6] eq "NA" and $info_id[9] eq "F"){

					$t5a3 = $t5a3 + $hash_isomirs{$key_isomir};

				}

			}

			# If it has 3' trimmings
			if ($info_id[3] ne "0") {

				# Search for the isomiRs that only have 3' trimmings 
				if ($info_id[2] eq "0" and $info_id[4] eq "0" and $info_id[5] eq "0" and $info_id[6] eq "NA" and $info_id[9] eq "F"){

					$t3 = $t3 + $hash_isomirs{$key_isomir};

				}

				# Search for the isomiRs that have both 3' trimmings and 5' additions
				elsif ($info_id[2] eq "0" and $info_id[4] ne "0" and $info_id[5] eq "0" and $info_id[6] eq "NA" and $info_id[9] eq "F"){

					$a5t3 = $a5a3 + $hash_isomirs{$key_isomir};

				}

				# Search for the isomiRs that have both 3' trimmings and 5' trimmings
				elsif ($info_id[2] eq "0" and $info_id[4] eq "0" and $info_id[5] ne "0" and $info_id[6] eq "NA" and $info_id[9] eq "F"){

					$t5t3 = $t5t3 + $hash_isomirs{$key_isomir};

				}

			}

			# Search for the isomiRs that only have 5' additions 
			elsif ($info_id[2] eq "0" and $info_id[3] eq "0" and $info_id[4] ne "0" and $info_id[5] eq "0" and $info_id[6] eq "NA" and $info_id[9] eq "F"){

				$a5 = $a5 + $hash_isomirs{$key_isomir};
			}

			# Search for the isomiRs that only have 5' trimmings 
			elsif ($info_id[2] eq "0" and $info_id[3] eq "0" and $info_id[4] eq "0" and $info_id[5] ne "0" and $info_id[6] eq "NA" and $info_id[9] eq "F"){

				$t5 = $t5 + $hash_isomirs{$key_isomir};
			}

			# Search for the isomiRs that have internal editings 
			elsif ($info_id[8] eq "Ed"){

				$editings = $editings + $hash_isomirs{$key_isomir};

			}

			# Search for the isomiRs that have tailings 
			elsif ($info_id[9] eq "T"){

				$tailings = $tailings + $hash_isomirs{$key_isomir};
			}
		}
	}

	
	# Array with the final counts to write in the output file
	my @counts = ($a5, $t5, $a3, $t3, $a5t3, $t5a3, $t5t3, $a5a3, $editings, $tailings);

	return @counts;

}

sub count_3p5p {

	# This subroutine receives a hash of counts and returns an array with the
	# the total counts for miRNAs that derive from the 3p end or 5p end 

	my (%hash_isomirs) = (@_);

	# Define counter variables
	my $p3 = 0;
	my $p5 = 0;

	foreach my $key_isomir (keys %hash_isomirs) {

		my @info_id = split("_", $key_isomir);
		my ($p3_p5) = $info_id[1] =~ /.*-(.*)/;


		if ($p3_p5 eq "3p") {

			$p3 = $p3 + $hash_isomirs{$key_isomir};
		}

		elsif ($p3_p5 eq "5p") {

			$p5 = $p5 + $hash_isomirs{$key_isomir};
		}

	}

	my @counts = ($p3, $p5);

	return @counts;

}


sub count_editings {

	my (%hash_isomirs) = (@_);
	
	my @typeOFeditings = ("A/U", "U/A", "A/G", "G/A", "A/C", "C/A", "G/C", "C/G");
	my %hash_editings;


	# Set the types of editings as keys and zero as default value
	for my $typeOf (@typeOFeditings) {
    	$hash_editings{$typeOf} //= 0;
	}

	foreach my $key_isomir (keys %hash_isomirs) {

		my @info_id = split("_", $key_isomir);
		my $editing = $info_id[6];
		my $Ed_SNP = $info_id[8];

		if ($Ed_SNP eq "Ed") {
			# For each type of editing in the editings hash,
			# if the found editing matches the editing on the editings hash
			# check if the value is zero. If so, change the value of the frequency
			# If the value != zero, append the value to the existing value
			# ex: if A/U =>0  : A/U => 2_3_1
			# ex if A/U => 1_5  : A/U => 1_5_2_3_1
			foreach my $type_editing (@typeOFeditings){
				if ($editing eq $type_editing) {
					if ($hash_editings{$editing} eq 0) {
						$hash_editings{$editing} = $hash_isomirs{$key_isomir};
					}
					else {
						$hash_editings{$editing} = $hash_editings{$editing} . '_' . $hash_isomirs{$key_isomir};
					}
				}
			}
		}
	}
	return %hash_editings;
	

}

sub count_tailings {

	my (%hash_isomirs) = (@_);

	my @typeOFtailings = ("UU", "UUU", "AA", "AAA", "GG", "GGG", "CC", "CCC");

	my %hash_tailings;


	# Set the types of tailings as keys and zero as default value
	for my $typeOf (@typeOFtailings) {
    	$hash_tailings{$typeOf} //= 0;
	}

	foreach my $key_isomir (keys %hash_isomirs) {

		my @info_id = split("_", $key_isomir);
		my $tailing = $info_id[10];

			# For each type of tailing in the tailings hash,
			# if the found tailing matches the tailing on the tailings hash
			# check if the value is zero. If so, change the value of the frequency
			# If the value != zero, append the value to the existing value
			# ex: if UU =>0  : UU => 2_3_1
			# ex if UU => 1_5  : UU => 1_5_2_3_1
			foreach my $type_tailing (@typeOFtailings){
				if ($tailing eq $type_tailing) {
					if ($hash_tailings{$tailing} eq 0) {
						$hash_tailings{$tailing} = $hash_isomirs{$key_isomir};
					}
					else {
						$hash_tailings{$tailing} = $hash_tailings{$tailing} . '_' . $hash_isomirs{$key_isomir};
					}
				}
			}
	
	}
	return %hash_tailings;

}


sub calculate_stDev {

	my (%hash_counts) = (@_);

	# Example of hash_counts:  miR-15b => 3_4_5_12_5

	# Calculates the standard deviation
	# std_dev = ( (sum (frequency - average) ** 2) / N) ** 0.5
	foreach my $hash_key (keys %hash_counts) {

		# Gets the sum of all the numbers in each value of the array
		my @numbers = split("_", $hash_counts{$hash_key});
		my $N = scalar @numbers;
		my $num_sum;

		foreach my $num (@numbers) {
			$num_sum += $num;
		}

		# Gets the average
		my $average = $num_sum / $N;
		my @number_for_sum;

		# Gets the square of each (mean_frequency - average)
		foreach my $value (@numbers) {
			my $value_for_sum = ($value - $average) ** 2;
			push(@number_for_sum, $value_for_sum);
		}

		# Gets the sum, to be divided by N
		my $value_sum;
		foreach my $val (@number_for_sum) {
			$value_sum += $val;
		}

		# Calculates the standard deviation 
		my $std_dev = ($value_sum / $N) ** 0.5;

		# Rounds the standard deviation to two decimal houses
		$std_dev = sprintf "%.2f", $std_dev;
		# Adds average_stDev to the value of each key
		$hash_counts{$hash_key} = $average . '_' . $std_dev;
	}

	return %hash_counts;
	# Ex: miR-15b => 5.8_3.19

}

sub get_scales {
	# Compare the files of two conditions to get the maximum value
	# (average + standard deviation) between both files, to set the
	# scale for the charts
	my (@file_array) = (@_);

	my $file_maxVal = 0;

	foreach my $file (@file_array) {

	# Open and read the files
	open (READ_FILE, "<$file") or die "Can't open $file: $!\n";

		# Go through each line and get the maximum value to set the scale
		#Example of a line: 15,157.00_87.00
		# We want the middle value to compare to $maxVal, in this case, 157.00
		while (defined (my $read_line = <READ_FILE>)) {
			chomp $read_line;

			my @lineElements = split(',', $read_line);	
				
			my @read_value = split('_', $lineElements[1]);
			# get the average + standard deviation
			my $read_avg_plus_stDev = $read_value[0] + $read_value[1];

			if ($read_avg_plus_stDev > $file_maxVal) {
				$file_maxVal = $read_avg_plus_stDev;
			}
		}
	}
	return sprintf("%.0f",$file_maxVal);
}

sub read_sizes {

	# This subroutine goes through an array of files,
	# opens each file and extracts the lenght ot the sequence
	# that aligned from the alignment file

	my (@all_files) = (@_);
	my $num_of_file = scalar @all_files;
	my %hash_sizes;

	foreach my $file (@all_files) {
		my %file_hash_sizes;
		open (GET_DATA_COUNTS, "<$file") or die "Can't open $file: $!\n";

		while (defined (my $line = <GET_DATA_COUNTS>)){

			chomp $line;

			# Skip header lines
			$line =~ /^@/ and next;
			my (@line_content) = split("\t", $line);
			my $aligned_or_not = $line_content[1];

			# Check if the read aligned
			if ($aligned_or_not eq "0" or $aligned_or_not eq "16") {

				# Extrat the number before the M character, that corresponds to the size of the read
				my $read_size = substr($line_content[5], 0, index($line_content[5], 'M'));

				if (exists $file_hash_sizes{$read_size}){

					$file_hash_sizes{$read_size} = $file_hash_sizes{$read_size} + 1;
				}

				else {

					$file_hash_sizes{$read_size} = 1

				}	
			}
		}

		# Append the values of each file to the final hash, separated by underscores
		# (to later calculate average and standard deviation)
		foreach my $key (keys %file_hash_sizes) {
			if (exists $hash_sizes{$key}) {
				$hash_sizes{$key} = $hash_sizes{$key} . '_' . $file_hash_sizes{$key};
			}
			else {
				$hash_sizes{$key} = $file_hash_sizes{$key};
			}
		}

	}

	return %hash_sizes;
}

sub calculate_positions {

	# Read an array of not_ambiguous files with isomiRs and
	# calculate the positions where the changes (regarding canonical)
	# occured. Produce a file for a chart

	# isomiR line denomination and example
	#Sequence-name-3pAddition/trimming_5pAddition/trimming_Editing/SNP/Alt(NucleotideChange_Position_TypeOf)_Tailing(T/F_TypeOf)
	#UCACAGUGAACCGGUCUCUUUAUC_ssc-miR-128_3_0_0_0_U/A_22_Alt_F_NA	3964

	# These types of changes will be taken into account:
	# Templated additions
	# Templated deletions
	# SNPs (A,U,G,C,ADAR[A/G]/PPR)
	# dbSNP (A,U,G,C)
	# Tailings


	my ($final_filename,$kingdom,@not_ambiguous_files) = (@_);

	# Get the number of files to later calculate the average of each count/position/typeOfChange
	my $number_of_files = scalar @not_ambiguous_files;

	my %SNP_positions;
	my %ED_positions;
	my %p5_aditions;
	my %p5_deletions;
	my %p3_aditions;
	my %p3_deletions;
	my %tailings;

	my $editing_type = "RNA Editing";

	#if ($kingdom eq 'animal') { $editing_type = "ADAR";}
	#elsif ($kingdom eq 'plant') { $editing_type = "PPR";}

	foreach my $isomir_file (@not_ambiguous_files) {
		# Open the not ambiguous isomiRs file
		open (GET_FILE_DATA_ISOMIRS,"<$isomir_file") or die "Can't open $isomir_file: $!\n";

		# Run through the not ambiguous isomiRs file
		while (defined (my $line = <GET_FILE_DATA_ISOMIRS>)) {
			chomp $line;

			# Split the line into separate values
			my @line_data = split("\t", $line);
			my $line_count = $line_data[1];
			my ($seq, $mir_name, $p3_adition, $p3_deletion, $p5_adition, $p5_deletion, $ed_nucleotide_pair, $ed_position, $ed_type, $tail_exists, $tail_type) = split('_', $line_data[0]);

			# Check to see if it not a line with multiple precursors (it would be in the place of the 3pA would be)
			if ($p3_adition eq 'multiple-precursors') {
				;
			}
			else {

				# The positions of some changes will have to be altered
				# Example of final positions:
				# -2,-1,5p,1,2,3,(...),21,3p,+1,+2,+3
				$p3_adition += 24;
				$p5_adition = 0 - $p5_adition;
				$p3_deletion = 24 - $p3_deletion;

				# Get the nucleotide that substituted the original
				# Ex: G/A => A
				my $ed_nucleotide; 
				if ($ed_nucleotide_pair eq "NA") {
					$ed_nucleotide = "NA";

				}
				else {
					my @ed_nucleotide_letters = split("/", $ed_nucleotide_pair);
					$ed_nucleotide = $ed_nucleotide_letters[1];					
				}

				# Push the values and counts of 5p and 3p aditions and deletions into hashes
				# If the value is zero, it will be deleted later
				# Ex: (%5p_aditions) : 2 => 123 | this means that there are 123 counts with two added nucleotides on the 5p end

				# Test if there is a key in the %p5_aditions hash
				# if there is, add the counts to the existing value
				# if not, create a key/value with the number of aditions and the corresponding count
				if (exists $p5_aditions{$p5_adition}) {
					$p5_aditions{$p5_adition} += $line_count;
				}
				else {
					$p5_aditions{$p5_adition} = $line_count;
				}


				# Same as before
				if (exists $p5_deletions{$p5_deletion}) {
					$p5_deletions{$p5_deletion} += $line_count;
				}
				else {
					$p5_deletions{$p5_deletion} = $line_count;
				}

				# Same as before
				if (exists $p3_aditions{$p3_adition}) {
					$p3_aditions{$p3_adition} += $line_count;
				}
				else {
					$p3_aditions{$p3_adition} = $line_count;
				}

				# Same as before
				if (exists $p3_deletions{$p3_deletion}) {
					$p3_deletions{$p3_deletion} += $line_count;
				}
				else {
					$p3_deletions{$p3_deletion} = $line_count;
				}

				# Check if there are tailings
				# If so, check the hash to see if there is a key/value pair corresponding
				# to that length. If so, add the counts to the value. If not, create the key/value pair
				my $tail_length;
				my $tail_position;

				if ($tail_exists eq "T") {
					$tail_length = length($tail_type);
					$tail_position = 24 + $tail_length;

					if (exists $tailings{$tail_position}) {
						$tailings{$tail_position} += $line_count;
					}
					else {
						$tailings{$tail_position} = $line_count;
					}
				}


				# Check if the editing is a SNP
				# If so, add the position and nucleotide's count to the hash's value, if it exists
				# If not, create a key/value pair
				if ($ed_type eq "SNP") {
					if (exists $SNP_positions{$ed_position}{$ed_nucleotide}) {
						$SNP_positions{$ed_position}{$ed_nucleotide} += $line_count;
					}
					else {
						$SNP_positions{$ed_position}{$ed_nucleotide} = $line_count;
					}

				}

				# Check if the editing is characterized as 'Ed'
				# If so, check if there are any ADAR editings
				# Add the position and nucleotide/typeOFEditing's count to the hash's value, if it exists
				# If not, create a key/value pair
				elsif ($ed_type eq "Ed") {

					if ($kingdom eq 'animal' && $ed_nucleotide_pair eq "A/G") {
						$ed_nucleotide = "ADAR";
					}
					elsif ($kingdom eq 'plant'){
						if ($ed_nucleotide_pair eq 'T/U' || $ed_nucleotide_pair eq 'U/C') {
							$ed_nucleotide = 'PPR';

						}
					}

					if (exists $ED_positions{$ed_position}{$ed_nucleotide}) {
						$ED_positions{$ed_position}{$ed_nucleotide} += $line_count;
					}
					else {
						$ED_positions{$ed_position}{$ed_nucleotide} = $line_count;
					}

				}
			}
			
		}
		}
		close(GET_FILE_DATA_ISOMIRS);


	my @SNP_changes = ("A", "U", "G", "C");
	my @ED_changes = ("A", "U", "G", "C", $editing_type);

	my $min_range = min(keys %p5_aditions);
	my $max_range = max(keys %p3_aditions);

	# Remove the 5p_aditions that are counted as zero
	# Aka where there were no 5p_aditions, the script counted
	if (exists $p5_aditions{"0"}) {
		delete $p5_aditions{"0"};
	}

	# Remove the 5p_deletions that are counted as zero
	# Aka where there were no 5p_deletions, the script counted
	if (exists $p5_deletions{"0"}) {
		delete $p5_deletions{"0"};
	}

	# Remove the 3p_aditions counted on position 24
	# which do not exist but were counted anyway
	if (exists $p3_aditions{"24"}) {
		delete $p3_aditions{"24"};
	}

	# Remove the 3p_deletions counted on position 24
	# which do not exist but were counted anyway
	if (exists $p3_deletions{"24"}) {
		delete $p3_deletions{"24"};
	}


	# Add alterations as zero were there were no alterations found
	# (so we can have the hash filled with every position, even if the count is zero)
	foreach my $i ($min_range..$max_range) {
		# Add the 3p aditions to the hash
		if (exists $p3_aditions{$i}) {
			;
		}
		else {
			$p3_aditions{$i} = 0;
		}

		# Add the 3p deletions to the hash
		if (exists $p3_deletions{$i}) {
			;
		}
		else {
			$p3_deletions{$i} = 0;
		}

		# Add the 5p additions to the hash
		if (exists $p5_aditions{$i}) {
			;
		}
		else {
			$p5_aditions{$i} = 0;
		}

		# Add the 5p deletions to the hash
		if (exists $p5_deletions{$i}) {
			;
		}
		else {
			$p5_deletions{$i} = 0;
		}

		# Add the SNPs to the hash
		if (exists $SNP_positions{$i}) {
			foreach my $snp_nucleotides (@SNP_changes) {
				if (exists $SNP_positions{$i}{$snp_nucleotides}) {
					;
				}
				else {
					$SNP_positions{$i}{$snp_nucleotides} = 0;
				}
			}
		}
		else {
			foreach my $snp_nucleotides (@SNP_changes) {
				$SNP_positions{$i}{$snp_nucleotides} = 0;
			}
		}

		# Add the Editings to the hash
		if (exists $ED_positions{$i}) {
			foreach my $ed_nucleotides (@ED_changes) {
				if (exists $ED_positions{$i}{$ed_nucleotides}) {
					;
				}
				else {
					$ED_positions{$i}{$ed_nucleotides} = 0;
				}
			}
		}
		else {
			foreach my $ed_nucleotides (@ED_changes) {
				$ED_positions{$i}{$ed_nucleotides} = 0;
			}
		}

		if (exists $tailings{$i}) {
			;
		}
		else {
			$tailings{$i} = 0;
		}

	}


	# Write to the final file

	# Open the not ambiguous isomiRs positions file (to make the charts)	 
	open(POSITIONS_FILE, '>', $final_filename);


	# Declare the final arrays to write in the file
	my (@p3A, @p3D, @p5A, @p5D, @SNP_A, @SNP_U, @SNP_G, @SNP_C, @Ed_A, @Ed_U, @Ed_G, @Ed_C, @Ed_special, @Tlng);

	# Calculate the max value for the chart scale
	my @maximum_values;

	# Calculate the average of the counts between all files 
	# and round to the lowest value (floor)

	#3p additions
	foreach my $p3Adds (sort {$a <=> $b} keys %p3_aditions){
		push(@p3A, floor($p3_aditions{$p3Adds}/$number_of_files));
		push(@maximum_values, floor($p3_aditions{$p3Adds}/$number_of_files));
	}

	# 3p deletions
	foreach my $p3Dels (sort {$a <=> $b} keys %p3_deletions){
		push(@p3D, floor($p3_deletions{$p3Dels}/$number_of_files));
		push(@maximum_values, floor($p3_deletions{$p3Dels}/$number_of_files));
	}

	#5p additions
	foreach my $p5Adds (sort {$a <=> $b} keys %p5_aditions){
		push(@p5A, floor($p5_aditions{$p5Adds}/$number_of_files));
		push(@maximum_values, floor($p5_aditions{$p5Adds}/$number_of_files));
	}

	#5p deletions
	foreach my $p5Dels (sort {$a <=> $b} keys %p5_deletions){
		push(@p5D, floor($p5_deletions{$p5Dels}/$number_of_files));
		push(@maximum_values, floor($p5_deletions{$p5Dels}/$number_of_files));
	}

	# SNPs
	foreach my $x_snp (sort {$a <=> $b} keys %SNP_positions){
		push(@SNP_A, floor($SNP_positions{$x_snp}{'A'}/$number_of_files));
		push(@SNP_U, floor($SNP_positions{$x_snp}{'U'}/$number_of_files));
		push(@SNP_G, floor($SNP_positions{$x_snp}{'G'}/$number_of_files));
		push(@SNP_C, floor($SNP_positions{$x_snp}{'C'}/$number_of_files));

		push(@maximum_values, floor($SNP_positions{$x_snp}{'A'}/$number_of_files),floor($SNP_positions{$x_snp}{'U'}/$number_of_files),
			floor($SNP_positions{$x_snp}{'G'}/$number_of_files),floor($SNP_positions{$x_snp}{'C'}/$number_of_files));
	}

	# Editings
	foreach my $x_ED (sort {$a <=> $b} keys %ED_positions){
		push(@Ed_A, floor($ED_positions{$x_ED}{'A'}/$number_of_files));
		push(@Ed_U, floor($ED_positions{$x_ED}{'U'}/$number_of_files));
		push(@Ed_G, floor($ED_positions{$x_ED}{'G'}/$number_of_files));
		push(@Ed_C, floor($ED_positions{$x_ED}{'C'}/$number_of_files));
		push(@Ed_special, floor($ED_positions{$x_ED}{$editing_type}/$number_of_files));

		push(@maximum_values, floor($ED_positions{$x_ED}{'A'}/$number_of_files),floor($ED_positions{$x_ED}{'U'}/$number_of_files),
			floor($ED_positions{$x_ED}{'G'}/$number_of_files),floor($ED_positions{$x_ED}{'C'}/$number_of_files),floor($ED_positions{$x_ED}{$editing_type}/$number_of_files));

	}


	# Tailings
	foreach my $tails (sort {$a <=> $b} keys %tailings){
		push(@Tlng, floor($tailings{$tails}/$number_of_files));
		push(@maximum_values, floor($tailings{$tails}/$number_of_files));
	}


	my @positions_to_write;
	# Change position 0 to 5p, 22 to 3p and add a plus sign to positions after 22
	foreach my $i ($min_range..$max_range) {

		if ($i > 24) {
			$i = '+' . ($i - 24);
		}
		push(@positions_to_write, $i);
	}

	my $file_header = '#' . join(',',@positions_to_write);

	my $to_write = $file_header . "\n3p_additions,".join('_', @p3A) . "\n3p_deletions,".join('_', @p3D) . "\n5p_additions,".join('_', @p5A) . "\n5p_deletions,".join('_', @p5D) . 
					"\nA(dbSNP),".join('_', @SNP_A) . "\nU(dbSNP),".join('_', @SNP_U) . "\nG(dbSNP),".join('_', @SNP_G) . "\nC(dbSNP),".join('_', @SNP_C) . 
					"\nA(SNP),".join('_', @Ed_A) . "\nU(SNP),".join('_', @Ed_U) . "\nG(SNP),".join('_', @Ed_G) . "\nC(SNP),".join('_', @Ed_C) . "\n$editing_type,".join('_', @Ed_special) . "\nTailings,".join('_',@Tlng);

	#print $to_write;

	# Write to file and close it
	print POSITIONS_FILE $to_write;
	close(POSITIONS_FILE);

	# Get the max val
	my $max_val = max(@maximum_values);
	return $max_val;
}




sub average {
	my @data = (@_);
	my $total=0;
	foreach my $value (@data){
		$total += $value;
	}
	my $average = $total/scalar @data;
	$average = sprintf "%.2f", $average;
	return $average;
}

sub stdev {
	my @data = (@_);
	if (scalar @data == 1) {
		return 0;
	}
	my $average = average(@data);
	my $sqtotal = 0;

	foreach my $value (@data) {
		$sqtotal += ($average-$value)**2;
	}
	my $stdev = ($sqtotal / (scalar @data-1)) ** 0.5;
	$stdev = sprintf "%.2f", $stdev;
	return $stdev;
}