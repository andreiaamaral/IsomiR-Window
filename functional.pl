#!/usr/bin/perl
use warnings;
use strict;

require "Species.pl";

# Receives the session's ID and the path to the input file
# The functional analysis may be integrated from the complete analysis or may start from an independent analysis
# When the file is provided by the user, it has to come in the FASTA format.
# When the file derives from the pipeline, it comes in the .TXT format

my($ID, $path_to_KB, $path_input_functional, $path_folder_inputs, $path_to_output, $species, $TargetScan, $miRanda, $targetFinder) = @ARGV;

my $path_input_topGO = "$path_folder_inputs"."/"."input_topgo.txt";
my $species_accession = Species::species_to_column($species, "genome_assembly");
# Obtain the species code
my $species_code = Species::species_to_column($species, "taxon_ID");


# If the user selects only TargetScan as target predictor
if (defined $TargetScan and $TargetScan != 0 and (!defined $miRanda or $miRanda == 0)) {# and (!defined $PITA  or $PITA == 0)){
	my $output_TargetScan = run_TargetScan($path_input_functional);
	system("cp $output_TargetScan $path_input_topGO");

}

# If the user selects only miRanda as target predictor
elsif (defined $miRanda and $TargetScan == 0 and $miRanda != 0){# and (!defined $PITA  or $PITA == 0)){

	my $output_miRanda = run_miRanda($path_input_functional);
	system("cp $output_miRanda $path_input_topGO");

}


# If the user selects TargetScan and miRanda as target predictors
elsif (defined $miRanda and $TargetScan != 0 and $miRanda != 0){# and (!defined $PITA  or $PITA == 0)){

	my $output_TargetScan = run_TargetScan($path_input_functional);
	my $output_miRanda = run_miRanda($path_input_functional);

	# Merge the outputs from the target predictors and selecting the unique and common genes
	system("sort $output_TargetScan $output_miRanda | uniq > $path_input_topGO");
}

# If the user selects TargetFinder (for plant species)
elsif ($targetFinder != 0) {
	my $output_TargetFinder = run_targetfinder($path_input_functional);
	system("cp $output_TargetFinder $path_input_topGO");
}

#####################################

# Path to R script to run topGO 
my $path_R_script = "R_topGO.R";


# Path to UTR universe
my $path_UTR_universe = $path_to_KB . "/".$species_accession."_UTR_universe.txt";
my $bioconductor_db = Species::species_to_column($species, "topGO_mapping");

my $no_topGO_file = "$path_to_output/no_TopGO.txt";

# Check to see if this species has a bioconductor db. If so, run TopGO
if ($bioconductor_db ne 'NONE') {
# Check if the file is empty
	if ( -z $path_input_topGO ) {
	    
	    # Define aqui o sitio para onde o queres mandar e o nome do ficheiro de erro
	    my $error_empty_file = "$path_to_output"."/"."ErrorFunctional.txt"; 

		# Create the input file for miRanda
		open(ERR_EMPTY, ">", "$error_empty_file") or die("Cannot write to file: $!\n");
		print ERR_EMPTY "The input file for TopGO is empty!\n";
	}

	# If it is not empty, run the R script
	else {

		system("Rscript --vanilla $path_R_script $ID $path_input_topGO $path_UTR_universe $path_to_output $bioconductor_db");

	}
}

# Write a file telling the user there is no file with the topGO results
else {
	open(NO_TOPGO, ">", "$no_topGO_file") or die("Cannot write to file: &!\n");
	print NO_TOPGO "There is no Bioconductor DB associated to this species, therefore there are no GO terms associated.\n";
	close(NO_TOPGO);

}

############################################### SUBROUTINES ####################################################



sub run_TargetScan {

	# Receives an input file
	my($input_file) = (@_);

	# Open the input file given in the argument
	open (GET_FILE_DATA_FUNC,"<$input_file") or die "Can't open $input_file: $!\n";

	# Declare the path of the input file for TargetScan
	my $input_targetscan = "$path_folder_inputs"."/"."input_targetScan.txt"; 

	# Create the input file for TargetScan
	open(INPUT_TARGETSCAN, ">", "$input_targetscan") or die("Cannot write to file: $!\n");


	my $flag = "TXT";
	my $first_line = <GET_FILE_DATA_FUNC>;

	# Check if the file contains a flag indicating that it is in the FASTA or TXT format
	if (index($first_line, $flag) != -1) {

		# Run through the input TXT file and prepare the file to be an input for 
		# the target predictors
		while (defined (my $line = <GET_FILE_DATA_FUNC>)) {
		
			# Clean the string out of whitespaces
			my @isomir = split /\s/ , $line;

			# Extract the ID of the isomiR/miRNA
			my $id_isomir = $isomir[0];

			# Extract the sequence of the isomiR/miRNA
			my $seq = substr($line, 0, index($line, '_'));

			#Extract the seed sequence 
			my $seed_seq = substr($seq, 1, 7);
			$seed_seq =~ s/T/U/g;

			#print $species_code;

			# Print to the file that will be the input for TargetScan
			print INPUT_TARGETSCAN "$id_isomir\t$seed_seq\t$species_code\n";

		}
	}

	else {

		# Run through the input FASTA file and prepare the file to be an input for 
		# TargetScan

		# Counter to go through every two lines
		my $counter = 0;
		my $stringToWrite = "";


		seek GET_FILE_DATA_FUNC, 0, 0;
		while (defined (my $line = <GET_FILE_DATA_FUNC>)) {

			chomp $line;

			# Skip lines that are only blank spaces
			$line =~ /^\s*$/ and next;

			if ($counter == 0) {

				$line =~ s/^.//;
				$line =~ s/\s//;

				$stringToWrite .= "$line"."_";

				$counter++;

			}	

			elsif ($counter == 1) {

				$line =~ s/\n//;
				my $seed_seq = substr($line, 1, 7);
				$seed_seq =~ s/T/U/g;

				$stringToWrite .= "$seed_seq";

				my $ID = substr($stringToWrite, 0, index($stringToWrite, '_'));
				my ($seq) = $stringToWrite =~ /_\s*(.*)\s*$/;

				

				print INPUT_TARGETSCAN "$ID\t$seq\t$species_code\n";

				# Reset the counter and the string 
				$counter = 0;
				$stringToWrite = "";
		
			}


		}
	}

	# Define the path to the UTR file for TargetScan
	my $path_to_UTR_targetScan = $path_to_KB . "/".$species_accession."_UTR_sequences.txt";

	# Define the path to the output file from TargetScan 
	my $path_to_output_targetScan = "$path_to_output"."/"."targetscan_output.txt"; 

	#Define the path to the TargetScan tool
	my $path_to_targetScan = "targetscan_70.pl";
	
	# Run TargetScan 
	system("echo 'yes' | perl $path_to_targetScan $input_targetscan $path_to_UTR_targetScan $path_to_output_targetScan");

	# Define the path to the list with unique genes file from Targetscan
	my $gene_list = "$path_to_output"."/"."gene_list_targetScan.txt";

	# Retrieve the genes of interest from the output
	system("cut -f1 $path_to_output_targetScan |sed 1d | sort | uniq | cut -f1 -d. > $gene_list");

	return "$gene_list";
}



sub run_miRanda {

	# Receives an input file
	my($input_file) = (@_);

	# Open the input file given in the argument
	open (GET_FILE_DATA_FUNC,"<$input_file") or die "Can't open $input_file: $!\n";


	my $flag = "TXT";
	my $first_line = <GET_FILE_DATA_FUNC>;
	my $input_miRanda; 

	# Check if the file contains a flag indicating that it is in the FASTA or TXT format
	# If it is in the TXT format
	if (index($first_line, $flag) != -1) {

		# Declare the path of the input file for miRanda
		$input_miRanda = "$path_folder_inputs"."/"."input_miranda.txt"; 

		# Create the input file for miRanda
		open(INPUT_MIRANDA, ">", "$input_miRanda") or die("Cannot write to file: $!\n");


		# Run through the input TXT file and prepare the file to be an input for 
		# the target predictors
		while (defined (my $line = <GET_FILE_DATA_FUNC>)) {
		
			# Clean the string out of whitespaces
			my @isomir = split /\s/ , $line;

			# Extract the ID of the isomiR/miRNA
			my $id_isomir = $isomir[0];

			# Extract the sequence of the isomiR/miRNA
			my $seq = substr($line, 0, index($line, '_'));

			# Sustitute all the Ts for Us
			$seq =~ s/T/U/g;

			# Print to the file that will be the input for miRanda
			print INPUT_MIRANDA ">"."$id_isomir\n";
			print INPUT_MIRANDA "$seq\n";

		}
	}

	else{

		$input_miRanda = "$path_folder_inputs"."/"."input_miranda.txt";

		system("cp $input_file $input_miRanda");
	}


	# Define the path to the UTR file for miRanda
	my $path_to_UTR_miRanda = $path_to_KB . "/".$species_accession."_UTR_sequences.fasta";

	# Define the path to the output file from miRAnda 
	my $path_to_output_miRanda = "$path_to_output"."/"."miRanda_output.txt"; 

	# Run miRanda 
	system("miranda $input_miRanda $path_to_UTR_miRanda -strict -out $path_to_output_miRanda");
	
	# Define the path to the list with unique genes file from Targetscan
	my $gene_list_miRanda = "$path_to_output"."/"."gene_list_miRanda.txt";

	# Create a file with the genes of interest
	system("cat $path_to_output_miRanda | grep '^>' | grep -v '>>' | awk '\$3 > 150 {print \$0}' | awk '\$4 < -7 {print \$0}' | cut -f2 | sort | uniq | cut -f1 -d. > $gene_list_miRanda");

	return "$gene_list_miRanda";
}

sub run_targetfinder {
	my($input_file) = (@_);
	# Open the input file given in the argument
	open (GET_FILE_DATA_FUNC,"<$input_file") or die "Can't open $input_file: $!\n";

	# Create the input file for TargetScan
	#open(OUTPUT_TARGETFINDER, ">", "$output_targefinder") or die("Cannot write to file: $!\n");

	my $flag = "TXT";
	my $first_line = <GET_FILE_DATA_FUNC>;
	my @targetFinder_seqs;

	# Check if the file contains a flag indicating that it is in the FASTA or TXT format
	if (index($first_line, $flag) != -1) {

		# Run through the input TXT file and prepare the file to be an input for 
		# the target predictors
		while (defined (my $line = <GET_FILE_DATA_FUNC>)) {
		
			# Clean the string out of whitespaces
			my @isomir = split /\s/ , $line;

			# Extract the sequence of the isomiR/miRNA
			my $seq = substr($line, 0, index($line, '_'));

			#Extract the seed sequence 
			my $seed_seq = substr($seq, 1, 7);
			$seed_seq =~ s/T/U/g;

			# Push to the array that will be the input for TargetScan;
			push @targetFinder_seqs, $seq;

		}
	}

	else {
		# Run through the input FASTA file and prepare the file to be an input for TargetFinder

		while (defined (my $line = <GET_FILE_DATA_FUNC>)) {
			chomp $line;
			# Skip lines that are only blank spaces
			$line =~ /^\s*$/ and next;
			my $first_char = substr $line,0;
			if ($first_char eq '>') {
				;
			}
			else {
				push @targetFinder_seqs, $line;
			}		
		}
	}

	my $path_to_UTR_targetFinder = $path_to_KB . "/".$species_accession."_UTR_sequences.fa";
	
	# Define the path to the output file from TargetScan 
	my $path_to_output_targetFinder = "$path_to_output"."/"."targetfinder_output.txt"; 

	#Define the path to the TargetScan tool
	my $path_to_targetFinder = "targetfinder.pl";

	# Run targetFinder
	foreach my $targetFinder_line (@targetFinder_seqs) {
		system("perl $path_to_targetFinder -s $targetFinder_line -d $path_to_UTR_targetFinder -r -p gff >> $path_to_output_targetFinder");
	}

	# Define the path to the list with unique genes file from TargeFinder
	my $gene_list = "$path_to_output"."/"."gene_list_targetFinder.txt";

	# Retrieve the genes of interest from the output
	system("cut -f1 $path_to_output_targetFinder |sort |uniq |grep -v 'No' > $gene_list");

	return "$gene_list";
}