
package Species;
use strict;
use warnings;


sub species_to_column {
	# Recieves a species scientific name and a column description and returns that column's value for the species
	# ex.: species_to_column("Sus_scrofa", "genome_assembly") returns "Sscrofa11.1"

	my ($species, $column) = (@_);

	my %column_number = (
		"miRDeep2_name" => 1,
		"genome_assembly" => 2,
		"miRBase_code" => 3,
		"taxon_ID" => 4,
		"topGO_mapping" => 5,
		"SNPs_database" => 6,
		"kingdom" => 7,
		);

	# Open the file with the Species information
	my $file = "species_hash.txt";
	open(INFO, '<', $file) or die "Could not open $file: $!";

	my $header_to_skip = <INFO>;  #skips the header

	my %hash_of_species;

	# Assign the species cientific name (row[0]) as a key to the column value
	foreach my $line (<INFO>) {
		chomp $line;
		my @row = split(";", $line);
		my $wanted_column = $column_number{$column};
		$hash_of_species{$row[0]} = $row[$wanted_column];

	}

	close(INFO);

	return $species = $hash_of_species{$species};

}

1;