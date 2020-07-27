#!/usr/bin/perl
use warnings;
use strict;
use 5.010;

require "./Species.pl";

my($path_to_KB, $path_to_condition, $path_bowtie_stats_file, $path_processing_files, $path_done_files, $path_filtered_SAM, $path_individual_results, $path_individual_ncrnas, $path_counts_ncrnas, $name_fastqfile, $number_mismatches, $number_multiple_hits, $specie) = @ARGV;


# Define the path of the input fastq file
my $fastqfile = "$path_to_condition"."/"."$name_fastqfile";
my $fastqfilelog = "$path_bowtie_stats_file"."/"."$name_fastqfile".".log";


# Check if the input fastq file is in the correct format
if (check_fastq_format($fastqfile) == 0){
	print "the fastq file is in the correct format\n";
}

# Check if the input fastq file has duplicate IDs, write the error
# # to a file and exit the script
elsif (check_fastq_format($fastqfile) == 1){

	my $error_output = "$path_processing_files"."/"."Error_format_IDs_$name_fastqfile.txt";
	open(FH, "> $error_output") or die ("Cannot write to file: $!\n");
	print FH "The file $name_fastqfile has duplicate IDs!\n" and die();
}

# Check if the input fastq file has not a suitable number of lines, 
# write the error to a file and exit the script
else {
	my $error_output = "$path_processing_files"."/"."Error_format_lines_$name_fastqfile.txt";
	open(FH, ">", "$error_output") or die("Cannot write to file: $!\n");
	print FH "The file $name_fastqfile has not a suitable number of lines (not multiple of 4)!\n" and die();
}


# Define the path of the output sam file for after running Bowtie
my $genome_assembly = Species::species_to_column($specie, "genome_assembly");
my $samfile = "$path_individual_results"."/"."$name_fastqfile.sam";
my $species_genome = $path_to_KB . "/" . $genome_assembly;


# Run the alignment and check if it ran

if (run_alignment("$number_mismatches", "$number_multiple_hits", "$species_genome", "$fastqfile", "$samfile") eq "0") {

	while (check_SAM_file($fastqfile, $samfile) != 0) {
		run_alignment("$number_mismatches", "$number_multiple_hits", "$species_genome", "$fastqfile", "$samfile");
	}

	my $Bowtie_finished_file = "$path_processing_files"."/"."Bowtie_finished_$name_fastqfile.txt";
	open(FH, ">", "$Bowtie_finished_file") or die("Cannot write to file: $!\n");
	print FH "The alignment with Bowtie is finished and the output SAM file is in the correct format!\n";
}
else {
	my $error_file = "$path_processing_files"."/"."Error_$name_fastqfile.txt";
	open(FH, ">", "$error_file") or die("Cannot write to file: $!\n");
	print FH "Bowtie could not run! Check if the files has no duplicate IDs\n" and die();
	
}


my $output_samfile_HTSeq = "$path_individual_results"."/"."HTSeq_$name_fastqfile.sam";


my $gff3_RNACentral = $path_to_KB . "/".$genome_assembly."_other_ncRNAs.gff3";
my $counts_HTSeq = "$path_individual_ncrnas"."/"."HTSeq_counts_$name_fastqfile.txt";


# Run HTSeq
run_HTSeq($samfile, $output_samfile_HTSeq, $gff3_RNACentral, $counts_HTSeq);


# Define the path and create the file with the general counts of the ncRNAs
my $general_counts_HTSeq_filename = "$path_counts_ncrnas"."/"."Counts_ncRNAs_HTSeq_$name_fastqfile.txt";
count_ncRNAs($counts_HTSeq, $general_counts_HTSeq_filename);


# Define the path of the SAM file with potential miRNAs that will be used in the following script
my $SAM_filtered_filename = "$path_filtered_SAM"."/"."Filtered_SAM_$name_fastqfile.sam";


if (extract_potential_miRNAs($output_samfile_HTSeq, $samfile, $SAM_filtered_filename) == 0) {
	
	my $done_file = "$path_done_files"."/"."Done_Sam_$name_fastqfile.txt";
	
	open(FH, ">", "$done_file") or die("Cannot write to file: $!\n");
	print FH "The filtered SAM file with potential miRNAs is done.\n";

}


else {
	
	my $error_output = "$path_processing_files"."/"."Error_no_SAMfile_$name_fastqfile.txt";
	open(FH, ">", "$error_output") or die("Cannot write to file: $!\n");
	print FH "The filtered SAM file for the $name_fastqfile was not created!\n" and die();
}

###################################### SUBROUTINES ######################################

sub check_fastq_format {

	my ($input_fastqfile) = (@_);

	#open the given file or show a message of fail attempt to open the specified file
	open (GET_FILE_DATA,"<$input_fastqfile") or die "Can't open $input_fastqfile: $!\n";

	my @unique;
	my %seen;
	my $number_of_lines = 0;

	while (defined (my $line = <GET_FILE_DATA>)){
		chomp $line;
		$number_of_lines++;

		if ((substr $line, 0) eq '@') {
			if (exists $seen{$line}) {
				return 1;
			}
			else {
				$seen{$line} = $line;
			}
		}
	}

	#close the file
	close GET_FILE_DATA;

	#check if the file has an appropriate number of lines (multiple of 4)
	if ($number_of_lines % 4 != 0){
		return 2;
	} 
	
	#@unique = ();
	#undef %seen;

	return 0;
}


sub run_alignment {

	my ($mismatches, $multiple_hits, $species, $input_fastqfile, $output_samfile) = (@_);

	if(system("bowtie -q -v$mismatches -m$multiple_hits -a --best --strata --sam $species $input_fastqfile $output_samfile 2>$fastqfilelog") == 0){

		return "0";
		
	}

}


sub check_SAM_file {
	# This subroutine checks if the output from bowtie (.sam) has the same number of 
	# sequences as the input fastq
	
	my ($input_fastqfile, $output_samfile) = (@_);

	use strict;
	use warnings;

	#open the given files or show a message of fail attempt to open the specified files
	open (GET_FILE_DATA_FASTQ,"<$input_fastqfile") or die "Can't open $input_fastqfile: $!\n";
	open (GET_FILE_DATA_SAM,"<$output_samfile") or die "Can't open $output_samfile: $!\n";

	my $number_of_IDs_fastq = 0;
	my $number_of_IDs_sam = 0;

	# Count the number of lines in fastq file
	my $countfq=0;
	$countfq++ while <GET_FILE_DATA_FASTQ>;
	$number_of_IDs_fastq=$countfq/4;

	# Go through every line of the SAM file
	my %seen= ();
	my @uniq=();
	my @list=();
	while (defined (my $line_sam=<GET_FILE_DATA_SAM>)){
		chomp $line_sam;
		my @fields=split ('\t', $line_sam);
		my $id_sam=$fields[0];
		if ($id_sam=~m/^@/){
		#	print "its header\n";
		}
		else {
			push (@list,$id_sam);
		}
	}
		

	foreach my $item (@list) {

	push(@uniq, $item) unless $seen{$item}++;
	}	
		 			 	
	$number_of_IDs_sam = scalar @uniq;
	
	#If the counts don't match, return 1
	if($number_of_IDs_fastq != $number_of_IDs_sam){
		return 1;
	}
	# If the counts match, return 0
	else {
		return 0;
	}
}


sub run_HTSeq {

	my ($input_samfile, $output_samfile, $input_gff3file, $output_counts_txt) = (@_);

	if (system("htseq-count --type=transcript --stranded=no --idattr=Name --samout=$output_samfile $input_samfile $input_gff3file > $output_counts_txt") == 0) {

		return "0";
	}
}

sub count_ncRNAs {

	#reads a txt file with the counts generated from HTseq
	my($counts_filename, $general_counts_filename) = (@_);

	use strict;
	use warnings;

	#open the given file or show a message of fail attempt to open the specified file
	open (GET_FILE_DATA,"<$counts_filename") or die "Can't open $counts_filename: $!\n";

	my $substring_notaligned = "not_aligned";
	my $substring_ambiguous = "ambiguous";

	my %count_hash;
	while (defined (my $line = <GET_FILE_DATA>)){
		chomp $line;

		#store the name and type in one variable, and the respective counts in another
		my ($name_type, $counts) = split("\t", $line);

		#if it is a line of a ncRNA
		if ($name_type =~ /^U/)  {
			my ($type) = $name_type =~ /=\s*(.*)\s*$/;
			$count_hash{$type}{total_ncRNAs} += $counts;
		}

		#if it is the line of not_aligned
		elsif (index($name_type, $substring_notaligned) != -1) {
			$count_hash{$name_type}{total_ncRNAs} += $counts;
		}

		#if it is the line of ambiguous
		elsif (index($name_type, $substring_ambiguous) != -1) {
			$count_hash{$name_type}{total_ncRNAs} += $counts;
		}
	}

	close(GET_FILE_DATA);
	
	#create a file with the counts for every type of ncRNA
	open(FH, "> $general_counts_filename") or die("Cannot write to file: $!\n");
	

	#write to the file
	foreach (sort keys(%count_hash)) {

		print FH "$_\t$count_hash{$_}{total_ncRNAs}\n";
	}

	return 0;
}

sub extract_potential_miRNAs {

	# Read a .sam file created after running HTSeq and the .sam used before running HTSeq
	my($sam_filename, $original_sam_filename, $filtered_SAM_filename) = (@_);

	use strict;
	use warnings;

	# Open the given sam files or show a message of fail attempt to open 
	#the specified files
	open (GET_FILE_DATA,"<$sam_filename") or die ("Can't open $sam_filename: $!\n");
	open (GET_FILE_DATA_O,"<$original_sam_filename") or die ("Can't open $original_sam_filename: $!\n");

	# Create the output file or show a message of fail attempt to 
	# create the specified file
	open(FH, "> $filtered_SAM_filename") or die ("Cannot write to file: $!\n");

	# Run through the first lines of the original .sam file and select the ones 
	# from the header
	
	while (defined (my $line = <GET_FILE_DATA_O> )){
		chomp $line;

		#write to the file if the line starts with a '@'
		if ($line =~ /^@/ ) {
			print FH "$line\n";
		}

		#stop running through the file if 
		else { 
			last;
		}	
	}

	# Run through every line of the .sam file from HTSeq and select all lines
	# except the __not_aligned, the ones only assigned to a rRNA and the
	# ambiguous only assigned to rRNAs
	
	while (defined (my $line2 = <GET_FILE_DATA> )){
		chomp $line2;
		my @cols = split("\t", $line2);
		# Get the string after the 'XF:Z:'
		# ex: from 'XF:Z:URS0000986EA5_type=rRNA' to 'URS0000986EA5_type=rRNA'
		my $alignmentFlag = substr $cols[-1], 5;

		# ex: URS00008B8881_type=rRNA
		my $RNAc_ID_start = substr $alignmentFlag, 0,3; # URS
		my $type_rRNA = substr $alignmentFlag, -10; # _type=rRNA

		my $ambiguous_Flag = 'False';

		my $rRNA_counter = 0;

		# Check if it is an ambiguous line. If so, check if it is only assigned fo rRNAs
		# example: XF:Z:__ambiguous[URS00008B8881_type=rRNA+URS0000237FE8_type=rRNA]
		# If this is the case, change the ambiguous_flag to True, as these lines are not to be considered
		if ((substr $alignmentFlag, 0,11) eq '__ambiguous'){

			# Get the types of RNA assigned to that line
			my @typesOfRNA = split /\+/,$alignmentFlag;

			foreach my $RNA (@typesOfRNA) {
				$RNA =~ s/]//; # remove the ] character from the end
				my $typeRNA = (split "=", $RNA)[1];
				if ($typeRNA eq 'rRNA') {
					$rRNA_counter++;
				}
			}
			# If the number of rRNAs is equal to the size of the array
			# aka if all the RNAs are rRNAs
			if (scalar @typesOfRNA == $rRNA_counter) {
				$ambiguous_Flag = 'True';
			}
		}

		# Check if it is a __not_aligned line or a line that is only assigned to a rRNA, or an ambiguous line only assigned to rRNAs
		# if so, do not print those lines
		# examples of (parts of) lines that are not to be printed:
		# XF:Z:__ambiguous[URS00008B8881_type=rRNA+URS0000237FE8_type=rRNA]
		# XF:Z:__not_aligned
		# XF:Z:URS00008B8881_type=rRNA

		if ($alignmentFlag eq '__not_aligned' || ($RNAc_ID_start eq 'URS' && $type_rRNA eq '_type=rRNA') || $ambiguous_Flag eq 'True') {
			; # pass and dont print those lines
		}

		# Examples of (parts of) lines that are going to be printed:
		# XF:Z:__ambiguous[URS00008B8881_type=rRNA+URS0000237FE8_type=lncRNA]
		# XF:Z:URS00008B8881_type=lncRNA

		else {
			splice @cols, 14, 1;
			print FH join("\t", @cols), "\n";
		}

	}

	close GET_FILE_DATA_O; 
	close GET_FILE_DATA;

	return 0;
}

1;