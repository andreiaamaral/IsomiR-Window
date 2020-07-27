#!/usr/bin/perl
use warnings;
use strict;
use Data::Dumper;
use Number::Range;
use Bio::Seq;
#use 5.010;

require "./Species.pl";

my($path_to_KB, $path_to_condition, $path_individual_results, $path_not_ambiguous, $path_ambiguous, $name_of_SAM_file, $species_genome) = @ARGV;


# Convert the species to the code of the reference genome
my $reference_genome_code = Species::species_to_column($species_genome, "genome_assembly");

my $genome_indexed_file = $path_to_KB . "/"."$reference_genome_code".".dna.fa";


# Create the index reference genome file if it has not yet been created
if (! -e $genome_indexed_file){
	system ("samtools faidx $genome_indexed_file");
}

# # Declare the path of the filtered SAM file that will be used
my $Filtered_SAM_name = "$path_to_condition"."/"."$name_of_SAM_file";

# Declare the path of the SNP_database that will be used
my $SNPs_file = Species::species_to_column($species_genome, "SNPs_database");
my $SNPs_database = $path_to_KB . "/"."$SNPs_file";


#print "$name_of_SAM_file\n";

my $output_samtools = "$path_individual_results"."/"."$name_of_SAM_file";

# Create the BAM file
system ("samtools view -b -S $Filtered_SAM_name > $output_samtools.bam");

system ("samtools sort -o $output_samtools.bam $output_samtools.bam");

#print "error";
# Declare the variables needed to execute the mpileup command
my $output_bam = "$output_samtools.bam";
my $mpileup_bcf = "$output_bam.bcf";
my $mpileup_final_bcf = "$output_bam"."_final".".bcf";
my $pathVCFutils = "vcfutils.pl";
my $SNPs_vcf = "$output_samtools"."_SNPs".".vcf";

# Run the mpileup command to create the vcf file with the SNPs 
run_mpileup($genome_indexed_file, $output_bam, $mpileup_bcf, $mpileup_final_bcf, $pathVCFutils, $SNPs_vcf);

# Remove the files that are not necessary
unlink glob "$path_individual_results/*.bcf";
unlink glob "$path_individual_results/*.csi";


# ################ Start classification of miRNAs #############################


#Use the created VCF to generate an annotated file through GATK

# Declare the path of the annotated SNPs file
my $SNPs_annotated = "$path_individual_results"."/".$name_of_SAM_file."_annotated_GATK.vcf";


# Run the GATK VariantAnnotator to annotate the vcf created from mpileup
system("java -jar GenomeAnalysisTK.jar -T VariantAnnotator -R $genome_indexed_file --variant $SNPs_vcf --dbsnp $SNPs_database  --alwaysAppendDbsnpId -o $SNPs_annotated");


# Open the annotated SNPs file
open (GET_FILE_DATA_SNP,"<$SNPs_annotated") or die "Can't open $SNPs_annotated: $!\n";


# Declare the hash to store the information needed from the annotated SNPs file
my %SNPsdata;
my @SNPs_data = ( );

# Run through the annotated SNPs file
while (defined (my $line4 = <GET_FILE_DATA_SNP>)) {

	# Skip header lines 
	$line4 =~ /^#/ and next;

	# Split the line into separate values
	my (@SNPs_data) = split("\t", $line4);

	# The key of the hash goes like: Chromossome_Position_Reference_Alteration -> id (. or rs..)
	$SNPsdata{"$SNPs_data[0]"."_"."$SNPs_data[1]"."_"."$SNPs_data[3]"."_"."$SNPs_data[4]"} = "$SNPs_data[2]"; 
}

#print Dumper(\%SNPsdata);

# Declare the path of the GFF file
my $gff3_miRBase =  $path_to_KB . "/"."$reference_genome_code"."_miRNAs.gff3";


# Open the GFF file
open (GET_FILE_DATA_GFF,"<$gff3_miRBase") or die "Can't open $gff3_miRBase: $!\n";


############# Create my database with the info from the GFF (miRBase) ##########
my %miRBasedata;
my @matureData =();

# Run through the gff file 
# Run through the gff file 
while (defined (my $line1 = <GET_FILE_DATA_GFF>)) {

	# Skip header lines and replace "\t" for ";"
    $line1 =~ /^#/ and next;
	$line1 =~ s/\t/;/g;
    chomp $line1;

    # Declare the variable that will differentiate if a line corresponds to a mature miRNA or primary transcript
    my $substr = "Derives_from";

    # Check if it is a line corresponding to a primary transcript miRNA and store its location as key and id as value into the hash
    if (index($line1, $substr) == -1) {
    	my ($chromossome, $source, $type, $start, $end, $score, $strand, $phase, $id, $alias, $name) = split(";", $line1);
    	my ($extracted_ID) = $id =~ /=\s*(.*)\s*$/; 
    	$miRBasedata{"$start"."_"."$end"."_"."$chromossome".";"."$strand"} = $extracted_ID;
			
	}

	# Check if it is a line corresponding to a mature miRNA and store
	# its location, primary transcript which it derives from, and name into an array
	else {
		my ($chromossome, $source, $type, $start, $end, $score, $strand, $phase, $id, $alias, $name, $primary_transcript) = split(";", $line1);
		my ($extracted_primary_transcript) = $primary_transcript =~ /=\s*(.*)\s*$/; 
		@matureData =($start, $end, $extracted_primary_transcript, $name);

		# Run through the keys of the hash and check if the id matches the primary transcript, push the mature data into the hash and delete the id
		foreach my $key (keys %miRBasedata){
			if ($matureData[2] eq ($miRBasedata{$key})){
				delete($miRBasedata{$key});
				$miRBasedata{$key} = [@matureData];
			}

			else{
				if ($matureData[2] eq $miRBasedata{$key}[2]){
					push @{ $miRBasedata{$key } }, @matureData; 
				}
				
				else {
					next;
				}	
			}
			
		}
		
	}

}

#print Dumper(\%miRBasedata);


# Open the SAM file
open (GET_FILE_DATA_SAM,"<$Filtered_SAM_name") or die "Can't open $Filtered_SAM_name: $!\n";


# Declare the hash and array that will store the info from the SAM file
# Declare the array that will store the reads that did not match to any miRNA
my %samData;
my $string_to_count = ();
my %remain;



my %frequencies;
while (defined (my $line2 = <GET_FILE_DATA_SAM>)) {
	$line2 =~ /^@/ and next;
	chomp $line2;

	my ($qname, $flag, $rname, $start_pos, $mapq, $cigar, $rnext, $pnext, $tlen, $seq, $qual, $XA, $MD, $NM) = split("\t", $line2);

	# Extract the lenght of the sequence and the end position of the sequence
	$cigar =~ tr/M//d;

	my $end_pos = 0;

	if ($flag == 0){

		$end_pos = $start_pos + $cigar -1;
	}

	elsif ($flag == 16) {

		$end_pos = $start_pos + $cigar -1;
	}
		

	# Extract the number of mismatches 
	my $num_mismatches = substr($NM, -1);

	# Extract the position of the mismatches
	my ($mismatches) = $MD =~ /Z:\s*(.*)\s*$/;

	# Define the string to count the frequency of the reads
	$string_to_count = "$rname,$mismatches,$start_pos,$end_pos,$seq,$num_mismatches,$flag";


	$frequencies{$string_to_count}++;

	
}


#print Dumper(\%frequencies);



foreach my $key (keys %frequencies) {

	my (@SAMdata) = split(",", $key);

	my $freq = "$frequencies{$key}";
	
	push (@SAMdata, $freq);

	$samData{"$SAMdata[2]"."_"."$SAMdata[3]"."_"."$SAMdata[4]".";"."$SAMdata[0]".";"."$SAMdata[6]"} = [@SAMdata];

}


#print Dumper(\%samData);



# Declare the path of the output files
my $output_not_ambiguous = "$path_not_ambiguous"."/"."$name_of_SAM_file"."_not_ambiguous.txt";
my $output_ambiguous = "$path_ambiguous"."/"."$name_of_SAM_file"."_ambiguous.txt";
my $handle_notAmb;

my $path_file_error = "$path_to_condition"."/"."Error_not_ambiguous.txt"; ###MUDAR AQUI O CAMINHO

if (!(open ($handle_notAmb,">$output_not_ambiguous"))) {

	open (OUTERROR,">$path_file_error");
	print OUTERROR "Can't open $output_not_ambiguous: $!\n";
}
#open ($handle_notAmb,">$output_not_ambiguous") or die "Can't open $output_not_ambiguous: $!\n";
open (OUTAMB,">$output_ambiguous") or die "Can't open $output_ambiguous: $!\n";


# Go through the SAM hash
foreach my $start_end_positions_SAM (keys %samData){

	# Set the variable to count the number of times the read matches a miRNA
	my $count=0;

	# Declare a hash which will contain all the possible matches for a sequence
	my %check_ambiguous;

	# Extract the start and end positions from the key
	my $start_in_key_SAM = substr($start_end_positions_SAM, 0, index($start_end_positions_SAM, '_'));
	my ($end_in_key_SAM) = $start_end_positions_SAM =~ /_([^_]+)_/;
	#my ($end_chr_SAM) = $start_end_positions_SAM =~ /;\s*(.*)\s*$/;
	my ($end_chr_SAM) = $start_end_positions_SAM =~ /;([^;]+);/;
	my $strand_SAM_16or0 = (split /;/, $start_end_positions_SAM)[-1];
	my $strand_SAM;
	my $chr_SAM = (split /;/, $start_end_positions_SAM)[-2];


	# Get the - or + information regarding the strand
	if ($strand_SAM_16or0 eq 16){
		$strand_SAM = '-';
	}
	elsif ($strand_SAM_16or0 eq 0) {
		$strand_SAM = '+';
	}

	# Go through the GFF hash
	foreach my $start_end_positions_miRBase (keys %miRBasedata){

		# Extract the start and end positions from the key
		my $start_in_key_miRBase = substr($start_end_positions_miRBase, 0, index($start_end_positions_miRBase, '_'));
		my ($end_in_key_miRBase) = $start_end_positions_miRBase =~ /_([^_]+)_/;
		my $strand_miRBase = (split /;/, $start_end_positions_miRBase)[-1];
		my $chr_miRBase = (split /_/, (split /;/, $start_end_positions_miRBase)[0])[-1];

		# Check if the read is contained in any pre miR
		if ($chr_SAM eq $chr_miRBase and $start_in_key_SAM >= $start_in_key_miRBase and $end_in_key_SAM <= $end_in_key_miRBase and $strand_SAM eq $strand_miRBase) {

			$check_ambiguous{$start_end_positions_miRBase} = $miRBasedata{$start_end_positions_miRBase};

			$count ++;
		}

	}

	# If it doesn't match anything, goes to the miRDeep2 Basket hash
	if ($count==0){

		$remain{$start_end_positions_SAM} = $samData{$start_end_positions_SAM};
	}

	# Get the size of the hash created for each element of the SAM 
	my $size_of_ash = keys %check_ambiguous;

	# It it has only one match (not ambiguous)
	if ($size_of_ash == 1) {

		# Select the only key from the hash
		foreach my $pre_mir( keys %check_ambiguous ){

			my $seq_miR = $samData{$start_end_positions_SAM}[4];
			my $flag_16_0 = $samData{$start_end_positions_SAM}[6];
			my $freq_miR = $samData{$start_end_positions_SAM}[7];
			my ($strand_premiR) = $pre_mir =~ /;\s*(.*)\s*$/;

			# Define the limit of overlap to consider a match
			#my $overlap_limit = (length($seq_miR))/2;


			# Check if it has 0 mismatches
			if ($samData{$start_end_positions_SAM}[5] == 0) {

				# No changes in nucleotides
				my $change_nucleotide = "NA"; #the change of the nucleotide (ex. T/A) or NA
				my $pos_mutation = 0; #0 or  the number of the position
				my $ed_or_SNP = "NA"; #editing, SNP or NA
				my $tailings = "F"; #true or false - T/F
				my $nucleotides_tail = "NA"; #the tailing (ex. AAA) or NA


				# Check if exists an overlap with the mature
				#if (overlap($check_ambiguous{$pre_mir}[0], $check_ambiguous{$pre_mir}[1], $start_in_key_SAM, $end_in_key_SAM) >= $overlap_limit){
				if (overlap($check_ambiguous{$pre_mir}[0], $check_ambiguous{$pre_mir}[1], $start_in_key_SAM, $end_in_key_SAM) > 0){

					
					my $difference_start = ($start_in_key_SAM - $check_ambiguous{$pre_mir}[0]);
					my $difference_end = ($end_in_key_SAM - $check_ambiguous{$pre_mir}[1]); 
					my $miRNA_name = $check_ambiguous{$pre_mir}[3];
					
					classification_isomiRs($handle_notAmb, $difference_start, $difference_end, $seq_miR, $miRNA_name, $freq_miR, $change_nucleotide, $pos_mutation, $ed_or_SNP, $tailings, $nucleotides_tail, $flag_16_0, $strand_premiR);

				}
				
				# If it does not match, check if the array contains the coordinates for the other mature miRNA
				elsif (scalar( @{ $check_ambiguous{$pre_mir} } ) == 8) {

					#Check if exists an overlap with the other mature miRNA
					#if (overlap($check_ambiguous{$pre_mir}[4], $check_ambiguous{$pre_mir}[5], $start_in_key_SAM, $end_in_key_SAM) >= $overlap_limit){
					if (overlap($check_ambiguous{$pre_mir}[4], $check_ambiguous{$pre_mir}[5], $start_in_key_SAM, $end_in_key_SAM) > 0){

						#my $seq_miR = $samData{$start_end_positions_SAM}[4];
						#my $freq_miR = $samData{$start_end_positions_SAM}[7];
						my $difference_start = ($start_in_key_SAM - $check_ambiguous{$pre_mir}[4]);
						my $difference_end = ($end_in_key_SAM - $check_ambiguous{$pre_mir}[5]); 
						my $miRNA_name = $check_ambiguous{$pre_mir}[7];


						classification_isomiRs($handle_notAmb, $difference_start, $difference_end, $seq_miR, $miRNA_name, $freq_miR, $change_nucleotide, $pos_mutation, $ed_or_SNP, $tailings, $nucleotides_tail, $flag_16_0, $strand_premiR)
						

					}

				} 

		
			}

			# Check if it has 1 mismatch
			elsif ($samData{$start_end_positions_SAM}[5] == 1) {


				# If it is a change in only one nucleotide, 
				# it doens't have tailings
				my $tailings = "F"; #true or false - T/F
				my $nucleotides_tail = "NA"; #the tailing (ex. AAA) or NA

				# Go to the MD code and extract the number before the letter to add to the 
				# start position and discover the exact location of the mutation 
				my $MD_code = "$samData{$start_end_positions_SAM}[1]";
				my ($nucleotides_to_add) = $MD_code =~ m/(\d+)[^\d*]/;
				my $position_mutation_sequence = $nucleotides_to_add + 1;

				# Add the extracted number to the start location
				my $pos_mutation = $start_in_key_SAM + $nucleotides_to_add + 2;

				# Extract the nucleotide in the MD code (nucleotide in reference)
				my ($nuc_reference) = $MD_code =~ m/\d+([^\d*]+)\d+/;

				if ($strand_premiR eq "-"){

					$nucleotides_to_add = $nucleotides_to_add + 1;
					$nuc_reference = rev($nuc_reference);
				}

				# Extrat the nucleotide in the read 
				#my $seq_miR = $samData{$start_end_positions_SAM}[4];
				my $nuc_read = substr($seq_miR, $nucleotides_to_add, 1);

				# Define the mutation in the read
				my $change_nucleotide = "$nuc_reference"."/"."$nuc_read";

				# Define a flag to know when to write
				my $found = "True";



				# Check if exists an overlap with the mature
				#if (overlap($check_ambiguous{$pre_mir}[0], $check_ambiguous{$pre_mir}[1], $start_in_key_SAM, $end_in_key_SAM) >= $overlap_limit){
				if (overlap($check_ambiguous{$pre_mir}[0], $check_ambiguous{$pre_mir}[1], $start_in_key_SAM, $end_in_key_SAM) > 0){
							
					my $difference_start = ($start_in_key_SAM - $check_ambiguous{$pre_mir}[0]);
					my $difference_end = ($end_in_key_SAM - $check_ambiguous{$pre_mir}[1]); 
					my $freq_miR = $samData{$start_end_positions_SAM}[7];
					my $miRNA_name = $check_ambiguous{$pre_mir}[3];
					my $lenght_sequence_miRBase = $check_ambiguous{$pre_mir}[1] - $check_ambiguous{$pre_mir}[0];

					my $check_mutation_flag = check_pos_mutation($seq_miR, $lenght_sequence_miRBase, $position_mutation_sequence, $strand_premiR);

					foreach my $SNP_location (keys %SNPsdata) {

						my @SNP_DATA = split /_/, $SNP_location;
						my $chromossome_SNP = $SNP_DATA[0];
						my $allele = $SNP_DATA[3];
						my $position_compare = $SNP_DATA[1];
						my $position_chrm_location = "$SNP_DATA[0]".":"."$SNP_DATA[1]";


						#print "$end_chr_SAM\t$chromossome_SNP\t$pos_mutation\t$position_compare\t$allele\t$nuc_read\n";

						if ("$end_chr_SAM" eq "$chromossome_SNP" and "$pos_mutation" eq "$position_compare" and "$allele" eq "$nuc_read"){


							if ($SNPsdata{$SNP_location} ne "."){

								my $ed_or_SNP = "SNP";

								if ($check_mutation_flag eq "False") {

									$ed_or_SNP = "NA";
									$position_mutation_sequence = 0;
									$change_nucleotide = "NA";
									
								}

								classification_isomiRs($handle_notAmb, $difference_start, $difference_end, $seq_miR, $miRNA_name, $freq_miR, $change_nucleotide, $position_mutation_sequence, $ed_or_SNP, $tailings, $nucleotides_tail, $flag_16_0, $strand_premiR);
								$found = "False";
							}

					
							else { 

								my $ed_or_SNP = "Ed";

								if ($check_mutation_flag eq "False") {

									$ed_or_SNP = "NA";
									$position_mutation_sequence = 0;
									$change_nucleotide = "NA";		
								}

								
								classification_isomiRs($handle_notAmb, $difference_start, $difference_end, $seq_miR, $miRNA_name, $freq_miR, $change_nucleotide, $position_mutation_sequence, $ed_or_SNP, $tailings, $nucleotides_tail, $flag_16_0, $strand_premiR);
								$found = "False";
							}
						}
						
					}

					if ($found eq "True" and $freq_miR > 1000) {

						my $ed_or_SNP =  "Alt";

						if ($check_mutation_flag eq "False") {

							$ed_or_SNP = "NA";
							$position_mutation_sequence = 0;
							$change_nucleotide = "NA";
						}

						classification_isomiRs($handle_notAmb, $difference_start, $difference_end, $seq_miR, $miRNA_name, $freq_miR, $change_nucleotide, $position_mutation_sequence, $ed_or_SNP, $tailings, $nucleotides_tail, $flag_16_0, $strand_premiR);

					}

				}
					

				# If it does not match, check if the array contains the coordinates for the other mature miRNA
				elsif (scalar( @{ $check_ambiguous{$pre_mir} } ) == 8) {

					#Check if exists an overlap with the other mature miRNA
					#if (overlap($check_ambiguous{$pre_mir}[4], $check_ambiguous{$pre_mir}[5], $start_in_key_SAM, $end_in_key_SAM) >= $overlap_limit){
					if (overlap($check_ambiguous{$pre_mir}[4], $check_ambiguous{$pre_mir}[5], $start_in_key_SAM, $end_in_key_SAM) > 0){

						my $difference_start = ($start_in_key_SAM - $check_ambiguous{$pre_mir}[4]);
						my $difference_end = ($end_in_key_SAM - $check_ambiguous{$pre_mir}[5]); 
						my $freq_miR = $samData{$start_end_positions_SAM}[7];
						my $miRNA_name = $check_ambiguous{$pre_mir}[7];
						my $lenght_sequence_miRBase = $check_ambiguous{$pre_mir}[5] - $check_ambiguous{$pre_mir}[4];

						my $check_mutation_flag = check_pos_mutation($seq_miR, $lenght_sequence_miRBase, $position_mutation_sequence, $strand_premiR);

						foreach my $SNP_location (keys %SNPsdata) {

							my @SNP_DATA = split /_/, $SNP_location;
							my $chromossome_SNP = $SNP_DATA[0];
							my $allele = $SNP_DATA[3];
							my $position_compare = $SNP_DATA[1];
							my $position_chrm_location = "$SNP_DATA[0]".":"."$SNP_DATA[1]";


							#print "$chromossome_SAM\t$chromossome_SNP\t$pos_mutation\t$location\t$allele\t$nuc_read\n";

							if ("$end_chr_SAM" eq "$chromossome_SNP" and "$pos_mutation" eq "$position_compare" and "$allele" eq "$nuc_read"){


								if ($SNPsdata{$SNP_location} ne "."){

									my $ed_or_SNP = "SNP";

									if ($check_mutation_flag eq "False") {

										$ed_or_SNP = "NA";
										$position_mutation_sequence = 0;
										$change_nucleotide = "NA";
										
									}

									classification_isomiRs($handle_notAmb, $difference_start, $difference_end, $seq_miR, $miRNA_name, $freq_miR, $change_nucleotide, $position_mutation_sequence, $ed_or_SNP, $tailings, $nucleotides_tail, $flag_16_0, $strand_premiR);
									$found = "False";
								}

					
								else { 

									my $ed_or_SNP = "Ed";

									if ($check_mutation_flag eq "False") {


										$ed_or_SNP = "NA";
										$position_mutation_sequence = 0;
										$change_nucleotide = "NA";
									}


									classification_isomiRs($handle_notAmb, $difference_start, $difference_end, $seq_miR, $miRNA_name, $freq_miR, $change_nucleotide, $position_mutation_sequence, $ed_or_SNP, $tailings, $nucleotides_tail, $flag_16_0, $strand_premiR);
									$found = "False";
								}
							}

						} 

						if ($found eq "True" and $freq_miR > 1000){

							my $ed_or_SNP =  "Alt";

							if ($check_mutation_flag eq "False") {

								$ed_or_SNP = "NA";
								$position_mutation_sequence = 0;
								$change_nucleotide = "NA";
							}

							classification_isomiRs($handle_notAmb, $difference_start, $difference_end, $seq_miR, $miRNA_name, $freq_miR, $change_nucleotide, $position_mutation_sequence, $ed_or_SNP, $tailings, $nucleotides_tail, $flag_16_0, $strand_premiR);
						}
					}
				}
			}


			# Check if it has 2 or 3 mismatches
			elsif ($samData{$start_end_positions_SAM}[5] == 2 or $samData{$start_end_positions_SAM}[5] == 3) {


				# Changes in more than one nucleotide can represent the 
				# existence of a tailing
				# Tailings only occur in the 3' end of the sequence
				# So, in order to be a tailing it has to be a repetition of 2 or 3
				# nucleotides in the end of the sequence, and never in the middle
				
				# Establish that the sequence doesn't have 
				# internal editings or SNPs
				my $MD_code = "$samData{$start_end_positions_SAM}[1]";
				my $change_nucleotide = "NA"; #the change of the nucleotide (ex. T/A) or NA
				my $pos_mutation = 0; #0 or the number of the position of the mutation
				my $ed_or_SNP = "NA"; #editing, SNP or NA
				

				# Check if exists an overlap with the mature
				#if (overlap($check_ambiguous{$pre_mir}[0], $check_ambiguous{$pre_mir}[1], $start_in_key_SAM, $end_in_key_SAM) >= $overlap_limit){
				if (overlap($check_ambiguous{$pre_mir}[0], $check_ambiguous{$pre_mir}[1], $start_in_key_SAM, $end_in_key_SAM) > 0){

					my $difference_start = ($start_in_key_SAM - $check_ambiguous{$pre_mir}[0]);
					my $difference_end = ($end_in_key_SAM - $check_ambiguous{$pre_mir}[1]); 
					my $miRNA_name = $check_ambiguous{$pre_mir}[3];
					my $tailings_or_not = check_tailings($seq_miR, $samData{$start_end_positions_SAM}[5], $MD_code, $flag_16_0, $difference_end, $difference_start);

					
					if ($tailings_or_not ne "NA"){


						my $tailings = substr($tailings_or_not, 0, index($tailings_or_not, '_'));
						my ($nucleotides_tail) = $tailings_or_not =~ /_\s*(.*)\s*$/;


						classification_isomiRs($handle_notAmb, $difference_start, $difference_end, $seq_miR, $miRNA_name, $freq_miR, $change_nucleotide, $pos_mutation, $ed_or_SNP, $tailings, $nucleotides_tail, $flag_16_0, $strand_premiR);
					}
					
				}

				
				# If it does not match, check if the array contains the coordinates for the other mature miRNA
				elsif (scalar( @{ $check_ambiguous{$pre_mir} } ) == 8) {

					#Check if exists an overlap with the other mature miRNA
					#if (overlap($check_ambiguous{$pre_mir}[4], $check_ambiguous{$pre_mir}[5], $start_in_key_SAM, $end_in_key_SAM) >= $overlap_limit){
					if (overlap($check_ambiguous{$pre_mir}[4], $check_ambiguous{$pre_mir}[5], $start_in_key_SAM, $end_in_key_SAM) > 0){

						#my $seq_miR = $samData{$start_end_positions_SAM}[4];
						#my $freq_miR = $samData{$start_end_positions_SAM}[7];
						my $difference_start = ($start_in_key_SAM - $check_ambiguous{$pre_mir}[4]);
						my $difference_end = ($end_in_key_SAM - $check_ambiguous{$pre_mir}[5]); 
						my $miRNA_name = $check_ambiguous{$pre_mir}[7];
						my $tailings_or_not = check_tailings($seq_miR, $samData{$start_end_positions_SAM}[5], $MD_code, $flag_16_0, $difference_end, $difference_start);
						


						if ($tailings_or_not ne "NA"){

							my $tailings = substr($tailings_or_not, 0, index($tailings_or_not, '_'));
							my ($nucleotides_tail) = $tailings_or_not =~ /_\s*(.*)\s*$/;

							classification_isomiRs($handle_notAmb, $difference_start, $difference_end, $seq_miR, $miRNA_name, $freq_miR, $change_nucleotide, $pos_mutation, $ed_or_SNP, $tailings, $nucleotides_tail, $flag_16_0, $strand_premiR);
						}
					}
				} 
			}
		}
	}

	elsif ($size_of_ash >= 2){

		my @matches = ( );
		my @mature_names = ( );
		
		foreach my $pre_mir ( keys %check_ambiguous ){

			push @matches, $check_ambiguous{$pre_mir}[2];
			my @mature_name = split(/_/, $check_ambiguous{$pre_mir}[3]);
			push @mature_names, (substr $mature_name[0], 5);
		}
		
		# Check if all the miRNAs in the mature_names list have the same name. If they have, it means that 
		# although they came from different precursors, they are the same mature miRNA.
		# If so, add a flag to the line
		
		# Check if the hash only has one key (miRNAs have the same name)
		my %string = map { $_, 1 } @mature_names;
		if (scalar (keys %string) == 1) {;
			my $mat_name = $mature_names[0];
			my $multiplePrecursor_IDs = join("_", @matches);

			# Sequence_matureName_'multiple-precursors'_matureID1_matureID2 frequency
			print $handle_notAmb $samData{$start_end_positions_SAM}[4] . "_" . $mat_name . "_multiple-precursors_" . $multiplePrecursor_IDs . "\t$samData{$start_end_positions_SAM}[7]\n";
		}
		
		else {

			# Sequence frequency IDs(pre-miRNAs)
			print OUTAMB "$samData{$start_end_positions_SAM}[4]\t$samData{$start_end_positions_SAM}[7]\t@matches\n";
		}

	}

}

#print Dumper(\%remain);

####################################### SUBROUTINES ###########################################


sub run_mpileup {

	my ($reference_genome_fasta, $output_filtered_sorted_miRNAs_bam, $output_filename_bcf, $output_final_bcf, $path_vcfutils, $SNPs_vcf) = (@_);

	system ("samtools-1.5/samtools mpileup --skip-indels -m 1 -E --BCF --output-tags DP,DV,DP4,SP -f $reference_genome_fasta -q20 -Q20 -o $output_filename_bcf $output_filtered_sorted_miRNAs_bam");
	system ("bcftools-1.5/bcftools index $output_filename_bcf indexed.bcf");
	system ("bcftools-1.5/bcftools call --skip-variants indels --multiallelic-caller --variants-only  --output-type b $output_filename_bcf -o $output_final_bcf");
	system ("bcftools-1.5/bcftools view $output_final_bcf | $path_vcfutils varFilter -a3 -d10 > $SNPs_vcf");

	return 0;
}

sub overlap {

	# This subroutine calculates the overlap between two sequences

	my ($start_position_miRBase, $end_position_miRBase, $start_position, $end_position) = (@_);

	my $range = Number::Range->new("$start_position_miRBase..$end_position_miRBase");
	my $range_to_test = Number::Range->new("$start_position..$end_position");
	my $overlap = 0;

	foreach my $int ($range_to_test->range) {
		
		if ( $range->inrange($int) ) {
			$overlap++;
		}

		else {
			next;
		}
	}
	return $overlap;
}


sub classification_isomiRs{

	my ($file, $start_sub, $end_sub, $sequence_sub, $miRNA_name_sub, $freq_sub, $nucleotide, $pos_chr, $ed_SNP, $tailing, $nucleotide_T, $flag, $strand_sub) = (@_);

	my ($miRNA_sub) = $miRNA_name_sub =~ /=\s*(.*)\s*$/;

	# Check if the position of the alteration 
	# is included in the sequence of the miRBase

	my $lenght_seq = length($sequence_sub);
	if ("$strand_sub" eq "-"){


		$sequence_sub = rev($sequence_sub);

		# Change the position of the mutation in the sequence, 
		# since it aligned in the other strand
		if ($pos_chr != 0) {   

			$lenght_seq = length($sequence_sub);
			$pos_chr = $lenght_seq - $pos_chr + 1;

		}

	}

	$nucleotide = T_to_U($nucleotide);
		
	$sequence_sub = T_to_U($sequence_sub); 

	# if the mutation is in the first or last nucleotide
	if ($pos_chr eq "1" or $pos_chr eq length($sequence_sub)) {
		;
	}

	else {
		# If both differences are positive, then we have an addition in 3' ends
		# and a trimming in 5' ends
		# Sequence_miRNA_3a_3t_5a_5t_(A/T)/NA_Pos_Ed/SNP_T/F_AAA/NA

		# Check the shift of the isomiR, relative to the canonical miRNA
		# Only consider isomiRs with <=2 nucleotides shift at 5p and <=5 shift nucleotides at 3p
		if (abs($start_sub) <= 2 and abs($end_sub) <= 5 ) {#and abs($start_sub) != 0 and abs($end_sub) != 0) {
			
			if ($start_sub == abs($start_sub) and $end_sub == abs($end_sub)) {

				if ("$strand_sub" eq "+") {

					print $file "$sequence_sub"."_"."$miRNA_sub".
					"_"."$end_sub"."_"."0"."_"."0"."_"."$start_sub"."_"."$nucleotide"."_"."$pos_chr"."_"."$ed_SNP"."_"."$tailing"."_"."$nucleotide_T\t$freq_sub\n";
				}

				elsif ("$strand_sub" eq "-"){

					print $file "$sequence_sub"."_"."$miRNA_sub".
					"_"."0"."_"."$start_sub"."_"."$end_sub"."_"."0"."_"."$nucleotide"."_"."$pos_chr"."_"."$ed_SNP"."_"."$tailing"."_"."$nucleotide_T\t$freq_sub\n";
				}
			}

			# If both differences are negative, then we have a trimming in 3' ends
			# and an addition in 5' ends
			elsif ($start_sub != abs($start_sub) and $end_sub != abs($end_sub)){

				$start_sub = abs($start_sub);
				$end_sub = abs($end_sub);

				if ("$strand_sub" eq "+") {

					print $file "$sequence_sub"."_"."$miRNA_sub".
					"_"."0"."_"."$end_sub"."_"."$start_sub"."_"."0"."_"."$nucleotide"."_"."$pos_chr"."_"."$ed_SNP"."_"."$tailing"."_"."$nucleotide_T\t$freq_sub\n";
				}

				elsif ("$strand_sub" eq "-"){

					print $file "$sequence_sub"."_"."$miRNA_sub".
					"_"."$start_sub"."_"."0"."_"."0"."_"."$end_sub"."_"."$nucleotide"."_"."$pos_chr"."_"."$ed_SNP"."_"."$tailing"."_"."$nucleotide_T\t$freq_sub\n";
				}

			}

			# If the difference in the start is positive, it is a trimming in 5' end
			# If the difference in the end is negative, it is also trimming in 3' end
			elsif ($start_sub == abs($start_sub) and $end_sub != abs($end_sub)){

				$end_sub = abs($end_sub);

				if ("$strand_sub" eq "+") {

					print $file "$sequence_sub"."_"."$miRNA_sub".
					"_"."0"."_"."$end_sub"."_"."0"."_"."$start_sub"."_"."$nucleotide"."_"."$pos_chr"."_"."$ed_SNP"."_"."$tailing"."_"."$nucleotide_T\t$freq_sub\n";
				}

				elsif ("$strand_sub" eq "-"){

					print $file "$sequence_sub"."_"."$miRNA_sub".
					"_"."0"."_"."$start_sub"."_"."0"."_"."$end_sub"."_"."$nucleotide"."_"."$pos_chr"."_"."$ed_SNP"."_"."$tailing"."_"."$nucleotide_T\t$freq_sub\n";
				}
			}

			# If the difference in the start is negative, it is an addition in 5' end
			# If the difference in the end is positive, it is also an addition in 3' end
			elsif ($start_sub != abs($start_sub) and $end_sub == abs($end_sub)){

				$start_sub = abs($start_sub);

				if ("$strand_sub" eq "+") {

					print $file "$sequence_sub"."_"."$miRNA_sub".
					"_"."$end_sub"."_"."0"."_"."$start_sub"."_"."0"."_"."$nucleotide"."_"."$pos_chr"."_"."$ed_SNP"."_"."$tailing"."_"."$nucleotide_T\t$freq_sub\n";
				}

				elsif ("$strand_sub" eq "-"){

					print $file "$sequence_sub"."_"."$miRNA_sub".
					"_"."$start_sub"."_"."0"."_"."$end_sub"."_"."0"."_"."$nucleotide"."_"."$pos_chr"."_"."$ed_SNP"."_"."$tailing"."_"."$nucleotide_T\t$freq_sub\n";
				}

		}}}

}

sub check_tailings{
	# This subroutine checks if the last 3 or 2 nucleotides can be considered a tailing
	my($sequence_test_tailings, $numberMM, $MDcode, $flag_16_0, $difference_end_sub, $difference_start_sub) = (@_);
	my $lenght_seq = length($sequence_test_tailings);
	my @all_nums = $MDcode =~ /(\d+)/g;

	if ($flag_16_0 eq "16"){

		$sequence_test_tailings = rev($sequence_test_tailings);
		$difference_end_sub = abs($difference_start_sub); 
	}


	if ($numberMM eq "2"){
		my $last_nucleotide1 = substr($sequence_test_tailings, -2, 1);
		my $last_nucleotide2 = substr($sequence_test_tailings, -1, 1);
		my $last_2_nucleotides = substr($sequence_test_tailings, -2);

		
		if ($last_nucleotide1 eq $last_nucleotide2 and $difference_end_sub >= 2){


			if ($all_nums[1] eq "0" and $all_nums[2] eq "0"){

				return "T"."_".T_to_U("$last_2_nucleotides");

			}
		}
	}

	elsif ($numberMM eq "3") {
		my $last_nucleotide1 = substr($sequence_test_tailings, -3, 1);
		my $last_nucleotide2 = substr($sequence_test_tailings, -2, 1);
		my $last_nucleotide3 = substr($sequence_test_tailings, -1, 1);
		my $last_3_nucleotides = substr($sequence_test_tailings, -3);
		my $last_2_nucleotides = substr($sequence_test_tailings, -2);

		if ($last_nucleotide1 eq $last_nucleotide2 and  $last_nucleotide2 eq $last_nucleotide3 and $difference_end_sub >= 3) {

			if ($all_nums[1] eq "0" and $all_nums[2] eq "0" and $all_nums[3] eq "0"){

				return "T"."_".T_to_U("$last_3_nucleotides");
			}
		}  

		elsif ($last_nucleotide2 eq $last_nucleotide3 and $difference_end_sub >= 2) {

			if ($all_nums[2] eq "0" and $all_nums[3] eq "0"){

				return "T"."_".T_to_U("$last_2_nucleotides");
			}

		}

	}

	return "NA";
}


sub rev{
	my ($read)=@_;
	my $revseq;
	#we take the sequence and make a new object Bio::Seq
		my $seqo= Bio::Seq->new( -seq => "$read");
		#print $seqo,"\n";
		#we reverse complement the sequence
		my $rev = $seqo->revcom();
		#we take out the sequence
		$revseq=$rev->seq;
		#print $revseq."\t".$read,"\n";
		return $revseq;		
	
}


sub T_to_U {

	# Receives a sequence with T's and turns it into U's
	my ($seq) = (@_);
	$seq =~ tr/T/U/; 
	return $seq;
}

sub change_code_alteration {

	# This subroutine gets as input a code like A/G and turns it into U/C

	my ($code_alteration) = (@_);

	my %hash_of_alterations = (
		"A/T"  => "U/A",
		"A/C"  => "U/G",
		"A/G"  => "U/C",
		"G/C"  => "C/G",
		"G/A"  => "C/U",
		"G/T"  => "C/A",
		"C/T"  => "G/A",
		"C/G"  => "G/C",
		"C/A"  => "G/U",
		"T/G"  => "A/C",
		"T/C"  => "A/G",
		"T/A"  => "A/U",
		
	);

	return $code_alteration = $hash_of_alterations{$code_alteration};
}

sub check_pos_mutation {

	# This sub checks if the mutation belongs to the sequence in the miRBase

	my ($seq_SAM, $seq_lenght_miRBase, $position_mutation, $strand) = (@_);

	my $True_or_False = "";

	if ($strand eq "+") {

		$position_mutation = $position_mutation - 2;
	}

	elsif ($strand eq "-") {

		my $seq_lenght_SAM = length($seq_SAM);
		$position_mutation = ($seq_lenght_SAM - $position_mutation) - 1;

	}

	if ($position_mutation <= $seq_lenght_miRBase) {

		$True_or_False = "True";

	}

	else {

		$True_or_False = "False";
	}

	return $True_or_False;

}