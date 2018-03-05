#!/usr/bin/perl --
#Copyright (c) 2014-2018, JF Flot <jflot@ulb.ac.be>
#
#Permission to use, copy, modify, and/or distribute this software for any purpose with or without fee is
#hereby granted, provided that the above copyright notice and this permission notice appear in all copies.
#
#THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES WITH REGARD TO THIS SOFTWARE
#INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
#ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS OF
#USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
#OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.

warnings;
use strict;
use Getopt::Std;
use List::Util 'shuffle';

my %opts;
getopts("i:n:r:s:vfgR", \%opts);
unless (defined($opts{i}) and (defined($opts{n}))) {print "Usage: perl countffrs.pl -i <input file in sequential FASTA format> -n <number of individuals to be sampled randomly from each simulated species> -s <number of simulated species; the program will divide the aligment into s equal blocks of sequences and consider each of them as coming from a different species, meaning that the haplotypes of heterozygotes will always be drawn from the same block> -r <number of replicate samplings> -v <verbiose> -f <save test dataset as FASTA> -g <save test dataset as gzipped FASTA> -R <generate Roehl output file for Network>\n"; exit};
unless (defined($opts{r})) {$opts{r}=1}; unless (defined($opts{s})) {$opts{s}=1};
#processing alignment
open(DATA, $opts{i}) || die "Error: $opts{i} cannot be opened!\n";
my@align=<DATA>;
    #remove all comments and endline characters
for (my $i=0; $i<=$#align; $i++)
	{chomp(@align);
	$/ = "\r";
	chomp(@align);
	$/ = "\n";
	if (substr($align[$i],0,1) eq ';') {splice(@align,$i, 1); $i=$i-1}};
    #makes a list of sequences
my @sequences;
for (my $i=0; $i<($#align+1)/2; $i++)
    {$sequences[$i]=$align[2*$i+1]};
    
    #repeat procedure several k times to make stats about the number of FFRs
for (my $k=0; $k<$opts{r}; $k++) {

	
	
	    #for each simulated species, pick up 2*$opts{n} sequences to form $opts{n} diploid individuals
    my @shuffled_indexes; my @picked_indexes; my @picked_sequences; my @table;
    for (my $s=0; $s<$opts{s}; $s++) { 
#print "$s\n";
    @shuffled_indexes=shuffle(($s*($#sequences+1)/$opts{s})..(((1+$s)*($#sequences+1)/$opts{s}-1))); 

#my $start=$s*($#sequences+1)/$opts{s}; my $end=(1+$s)*($#sequences+1)/$opts{s}-1; print "$start\t$end\n";
   
 	 @picked_indexes=@shuffled_indexes[0..2*$opts{n}-1];
 	 @picked_sequences=@sequences[@picked_indexes];

			     
 	for (my $i=0; $i<$opts{n}; $i++) {my $indiv = sprintf("%03d",$i+1); $table[($s*2*$opts{n})+2*$i][0]="T".($s+1)."indiv".$indiv."a"; $table[($s*2*$opts{n})+2*$i][1]=$picked_sequences[2*$i];
 							    $table[($s*2*$opts{n})+1+2*$i][0]="T".($s+1)."indiv".$indiv."b"; $table[($s*2*$opts{n})+1+2*$i][1]=$picked_sequences[1+2*$i]}}
 	
      
	
	
	    #for debugging
    #for (my $i=0; $i<($#table+1)/2; $i++) {print "$table[2*$i][0]\n$table[1+2*$i][0]\n"};
    
#write output dataset as FASTA if requested
	if ($opts{f} or $opts{g}) {
	    open (F, "> rep".$k."_".$opts{i}) || die "Cannot write in out.fasta!";
	    for (my $i=0; $i<($#table+1)/2 ; $i++) {print F ">".$table[2*$i][0]."\n".$table[2*$i][1]."\n>".$table[1+2*$i][0]."\n".$table[1+2*$i][1]."\n"};
	    close F};
    if ($opts{g}) {system("gzip -v9 rep".$k."_".$opts{i}) && die "unable to gzip output FASTA"} ;
	
	
	#all sequences are simulated so have supposedly the same length
    my $maxseqlength=length($table[0][1]);
	for (my$i=0; $i<=$#table; $i++)
		{my $l=length($table[$i][1]); unless ($l==$maxseqlength) {print "Not all sequences in the input dataset $opts{i} have the same length!"; exit}};
	
	#make hash of haplotypes
	my %haplotypes;
	for (my$i=0; $i<=$#table; $i++) {push(@{$haplotypes{$table[$i][1]}}, $table[$i][0])};
	my $number_of_haplotypes = keys %haplotypes;
	
	if ($opts{v}) {
	print "There are $number_of_haplotypes haplotypes in your dataset of $opts{n} individuals.\n";
	}
	
	#make table of haplotypes
	my @haplotable; my @haplotable2;
	my$count=0;
	foreach my$key (sort {$#{$haplotypes{$b}} <=>$#{$haplotypes{$a}}} keys %haplotypes) 
			{@{$haplotable[$count]}=@{$haplotypes{$key}}; $haplotable2[$count]=$key; $count++};
	
	#display list of haplotypes
	
	if ($opts{v}) {
	for (my$i=0; $i<$number_of_haplotypes; $i++)
			{print "Haplotype ", $i+1," is found ", $#{$haplotable[$i]}+1, " time(s) in the alignment, under the name(s): ", join(', ', @{$haplotable[$i]}), ".\n"};
	}
	
	#make heterozygotes table
	my @heterozygotes;
	for (my$i=0; $i<$#table; $i++)
	    {my $previous = $table[$i][0]; my $next = $table[$i+1][0]; my $endprevious=chop($previous); my 
	$endnext=chop($next);
	if (((($endprevious eq 'a') and ($endnext eq 'b')) or (($endprevious eq 'b') and ($endnext eq 'c')) or 
	(($endprevious eq 'c') and ($endnext eq 'd')) or (($endprevious eq 'd') and ($endnext eq 'e')) or 
	(($endprevious eq 'e') and ($endnext eq 'f'))) and (($previous eq $next) and ($table[$i][1] ne 
	$table[$i+1][1]))) {push(@heterozygotes, $previous) unless grep(/$previous/,@heterozygotes)}};
	
	#display list of heterozygotes
	
	if ($opts{v}) {
	print "A total of ", scalar(@heterozygotes), " individuals contain more than one haplotype: ", join(', ',@heterozygotes),".\n";
	}
	
	#make hash of connections
	my %connect;
	foreach my$indiv (@heterozygotes) 
		{my $ka; my $kb; my $kc; my $kd; my $ke; my $kf; for (my$i=0; $i<=$#haplotable; $i++)
			{if (grep(/$indiv+a/, @{$haplotable[$i]})) {$ka=$i+1}; if (grep(/$indiv+b/, 
	@{$haplotable[$i]})) {$kb=$i+1};
			 if (grep(/$indiv+c/, @{$haplotable[$i]})) {$kc=$i+1}; if (grep(/$indiv+d/, 
	@{$haplotable[$i]})) {$kd=$i+1}
			 if (grep(/$indiv+e/, @{$haplotable[$i]})) {$ke=$i+1}; if (grep(/$indiv+f/, 
	@{$haplotable[$i]})) {$kf=$i+1}};
	
			($ka,$kb,$kc,$kd,$ke,$kf)=sort {$a<=>$b} ($ka,$kb,$kc,$kd,$ke,$kf);
			$connect{"$ke and $kf"}{$indiv}='OK';
			if ($kd) {$connect{"$kd and $ke"}{$indiv}='OK'; $connect{"$kd and $kf"}{$indiv}='OK'};
			if ($kc) {$connect{"$kc and $kd"}{$indiv}='OK'; $connect{"$kc and $ke"}{$indiv}='OK'; 
	$connect{"$kc and $kf"}{$indiv}='OK'};
			if ($kb) {$connect{"$kb and $kc"}{$indiv}='OK'; $connect{"$kb and $kd"}{$indiv}='OK'; 
	$connect{"$kb and $ke"}{$indiv}='OK'; $connect{"$kb and $kf"}{$indiv}='OK'};
			if ($ka) {$connect{"$ka and $kb"}{$indiv}='OK'; $connect{"$ka and $kc"}{$indiv}='OK'; 
	$connect{"$ka and $kd"}{$indiv}='OK'; $connect{"$ka and $ke"}{$indiv}='OK'; $connect{"$ka and 
	$kf"}{$indiv}}};
	
	#remove spurious connections between one haplotype and itself (useful when some sites have been removed from the alignment, in which case individuals may have several times the same haplotype)
	for (my$i=0; $i<=$#haplotable; $i++) {delete($connect{"$i and $i"})};
	
	#display list of connections
	
	if ($opts{v}) {
	foreach my$key (sort {scalar(keys %{$connect{$b}}) <=> scalar(keys %{$connect{$a}})} keys %connect)
	{my @connectingindividuals = sort {$#{$connect{$a}} <=>$#{$connect{$b}}} keys %{$connect{$key}};
	print "Haplotypes $key co-occur in ", scalar(@connectingindividuals), " individual(s), namely: ", 
	join(', ', sort @connectingindividuals),".\n"};
	}
	
	#build FFRs
	foreach my$hetero(@heterozygotes) {my $keya; my $keyb; my $keyc; my $keyd; my $keye; my $keyf; my 
	@reunion;
	while (my ($k, $v) = each (%haplotypes)) {if (grep(/$hetero+a/, @{$v})) {$keya=$k}; if 
	(grep(/$hetero+b/, @{$v})) {$keyb=$k};
	if (grep(/$hetero+c/, @{$v})) {$keyc=$k}; if (grep(/$hetero+d/, @{$v})) {$keyd=$k}; if 
	(grep(/$hetero+e/, @{$v})) {$keye=$k};
	 if (grep(/$hetero+f/, @{$v})) {$keyf=$k}};
	if (($keyb) and ($keya ne $keyb)) {@reunion=(@{$haplotypes{$keya}},@{$haplotypes{$keyb}});
	@{$haplotypes{$keya}}=@reunion; delete($haplotypes{$keyb}) unless (($keyb eq $keyc) or ($keyb eq $keyd) 
	or ($keyb eq $keye) or ($keyb eq $keyf))};
	if (($keyc) and ($keya ne $keyc)) {@reunion=(@{$haplotypes{$keya}},@{$haplotypes{$keyc}});
	@{$haplotypes{$keya}}=@reunion; delete($haplotypes{$keyc}) unless (($keyc eq $keyd) or ($keyc eq $keye) 
	or ($keyc eq $keyf))}; 
	if (($keyd) and ($keya ne $keyd)) {@reunion=(@{$haplotypes{$keya}},@{$haplotypes{$keyd}});
	@{$haplotypes{$keya}}=@reunion; delete($haplotypes{$keyd}) unless (($keyd eq $keye) or ($keyd eq 
	$keyf))}; 
	if (($keye) and ($keya ne $keye)) {@reunion=(@{$haplotypes{$keya}},@{$haplotypes{$keye}});
	@{$haplotypes{$keya}}=@reunion; delete($haplotypes{$keye}) unless ($keye eq $keyf)};
	if (($keyf) and ($keya ne $keyf)) {@reunion=(@{$haplotypes{$keya}},@{$haplotypes{$keyf}});
	@{$haplotypes{$keya}}=@reunion; delete($haplotypes{$keyf})};
	};
	
	#remove duplicate sequence names in haplotypes hash
	foreach my $k (keys %haplotypes) {my %hash = map { $_, 1 } @{$haplotypes{$k}}; @{$haplotypes{$k}}= keys 
	%hash};
	
	
	#sort sequence names in FFRs
	foreach my$key (keys %haplotypes) {@{$haplotypes{$key}} = sort @{$haplotypes{$key}} };
	
	
	#display FFRs in terms of sequences
	
	if ($opts{v}) {
	if (scalar(keys %haplotypes) == 1) {print "There is a single FFR in your dataset.\n"}
	else {print "There are ", scalar(keys %haplotypes), " FFRs in your dataset.\n"}
	};
	
	if ($opts{v}) {
	foreach my$key (sort {$#{$haplotypes{$b}} <=> $#{$haplotypes{$a}}} keys %haplotypes) {print "One FFR contains the following ", scalar(@{$haplotypes{$key}})," sequence(s): ", join(', ',@{$haplotypes{$key}},".\n")}; 
	}
	my $total_number_sequences=0;
	foreach my$key (keys %haplotypes) {$total_number_sequences+= scalar(@{$haplotypes{$key}})};
	if ($opts{v}) {
	print "Total number of sequences in the dataset: $total_number_sequences.\n";
	}
	#display FFRs in terms of haplotypes
	if ($opts{v}) {print "Or, in terms of haplotypes:\n"};
	my %haplonumbers;
	for (my $i=0; $i<=$#haplotable; $i++) {foreach my $seq (@{$haplotable[$i]}) 
	{$haplonumbers{$seq}=$i+1}};
	
	my %haplotypes2;
	foreach my$key (keys %haplotypes) {foreach my $seq (@{$haplotypes{$key}}) 
	{$haplotypes2{$haplonumbers{$seq}}=$key}};
	
	my %haplotypes3;
	foreach my$key (keys %haplotypes2) {push(@{$haplotypes3{$haplotypes2{$key}}},,$key)};
	
	if ($opts{v}) {
	foreach my$key (sort {$#{$haplotypes3{$b}} <=>$#{$haplotypes3{$a}}} keys %haplotypes3) 
	{print "One FFR contains the following ", scalar(@{$haplotypes3{$key}})," haplotype(s): ", join(', ', 
	sort {$a <=> $b} @{$haplotypes3{$key}}),".\n"}; 
	}
	
	#display FFRs in terms of individuals	
	if ($opts{v}) {
	print "Or, in terms of individuals:\n";
	}
	
	foreach my$key (keys %haplotypes) {
	for (my$i=0; $i<$#{$haplotypes{$key}}; $i++)
	    {my $previous = ${$haplotypes{$key}}[$i]; my $next = ${$haplotypes{$key}}[$i+1]; my $nextnext = 
	${$haplotypes{$key}}[$i+2];
	my $nextnextnext = ${$haplotypes{$key}}[$i+3]; my $nextnextnextnext = ${$haplotypes{$key}}[$i+4]; 
	my $nextnextnextnextnext = ${$haplotypes{$key}}[$i+5];
	my $endprevious=chop($previous); my $endnext=chop($next); my $endnextnext=chop($nextnext); my 
	$endnextnextnext=chop($nextnextnext);
	my $endnextnextnextnext=chop($nextnextnextnext); my 
	$endnextnextnextnextnext=chop($nextnextnextnextnext);
	if (($endprevious eq 'a') and ($endnext eq 'b') and ($previous eq $next)) 
		{${$haplotypes{$key}}[$i]=$previous; splice (@{$haplotypes{$key}},$i+1,1);
		if (($endnextnext eq 'c') and ($previous eq $nextnext)) 
			{splice (@{$haplotypes{$key}},$i+1,1);
			if (($endnextnextnext eq 'd') and ($previous eq $nextnextnext)) 
				{splice (@{$haplotypes{$key}},$i+1,1);
				if (($endnextnextnextnext eq 'e') and ($previous eq $nextnextnextnext))
					{splice (@{$haplotypes{$key}},$i+1,1);
					if (($endnextnextnextnextnext eq 'e') and ($previous eq 
	$nextnextnextnextnext))
						{splice (@{$haplotypes{$key}},$i+1,1);
	}}}}}}};
	
	if ($opts{v}) {
	foreach my$key (sort {$#{$haplotypes{$b}} <=>$#{$haplotypes{$a}}} keys %haplotypes) {print "One FFR contains the following ", scalar(@{$haplotypes{$key}})," individual(s): ", (join(', ',@{$haplotypes{$key}})), "\n"};
	} 
    else { foreach my$key (sort {$#{$haplotypes{$b}} <=>$#{$haplotypes{$a}}} keys %haplotypes) {print join(' ',@{$haplotypes{$key}}).";"}; print "\n"};
	my $total_number_individuals=0;
	foreach my$key (keys %haplotypes) {$total_number_individuals+= scalar(@{$haplotypes{$key}})};
	
	if ($opts{v}) {
	print "Total number of individuals in the dataset: $total_number_individuals.\n\n";
	}
	
	if ($opts{R}) {
	#generate Roehl output
	if (!open (F, "> $opts{i}.rdf")){die "Writing error: $!"};
	
	my @haplotable3; my@varpos=();
	for (my$i=0; $i<=$#haplotable2; $i++) {@{$haplotable3[$i]}=split('', $haplotable2[$i])};
	
	POSITION: for (my$j=0; $j<=$maxseqlength; $j++) {
		for (my$i=0; $i<$#haplotable3; $i++) { if ($haplotable3[$i][$j] ne $haplotable3[$i+1][$j]) 
	{push(@varpos, $j+1); next POSITION}}};
	
	my@varposdigits; for (my$i=0; $i<=$#varpos; $i++) {@{$varposdigits[$i]}=split('',$varpos[$i])};
	for (my$i=0; $i<=$#varpos; $i++) {for (my $j=1+$#{$varposdigits[$i]}; $j<=$#{$varposdigits[$#varpos]}; 
	$j++) {$varposdigits[$i][$j]=' '}};
	print F '       '; for (my$i=0; $i<=$#varpos; $i++) {print F $varposdigits[$i][0]}; print F "\r\n";
	print F '       '; for (my$i=0; $i<=$#varpos; $i++) {print F $varposdigits[$i][1]}; print F "\r\n";
	print F '       '; for (my$i=0; $i<=$#varpos; $i++) {print F $varposdigits[$i][2]}; print F "\r\n";
	print F '       '; for (my$i=0; $i<=$#varpos; $i++) {print F $varposdigits[$i][3]}; print F "\r\n";
	print F '       '; for (my$i=0; $i<=$#varpos; $i++) {print F $varposdigits[$i][4]}; print F "\r\n";
	print F '       '; for (my$i=0; $i<=$#varpos; $i++) {print F $varposdigits[$i][5]}; print F "\r\n";
	
	my @haplonumbers; for (my$i=0; $i<$number_of_haplotypes; $i++) {$haplonumbers[$i]=1+$i};
	my @haplonumbersdigits;
	for (my$i=0; $i<=$#haplonumbers; $i++) {@{$haplonumbersdigits[$i]}=split('',$haplonumbers[$i])};
	for (my$i=0; $i<=$#haplonumbers; $i++) {for (my $j=1+$#{$haplonumbersdigits[$i]}; $j<=6; $j++) 
	{$haplonumbersdigits[$i][$j]=' '}};
	
	my @samplesizes; for (my$i=0; $i<=$#haplonumbers; $i++) {$samplesizes[$i]=$#{$haplotable[$i]}+1};
	my @samplesizesdigits;
	for (my$i=0; $i<=$#samplesizes; $i++) {@{$samplesizesdigits[$i]}=split('',$samplesizes[$i])};
	for (my$i=0; $i<=$#samplesizes; $i++) {if ($samplesizes[$i]>999) {print "One haplotype is present more than 999 times. As a result, no Network input file can be generated."; close F; exit} elsif ($samplesizes[$i]<10) {unshift(@{$samplesizesdigits[$i]},'  ')} elsif ($samplesizes[$i]<100) {unshift(@{$samplesizesdigits[$i]},' ')}};
	
	for (my$i=0; $i<=$#haplonumbers; $i++) {print F @{$haplonumbersdigits[$i]}; 
		for (my $j=0; $j<=$#varpos; $j++) {print F "$haplotable3[$i][$varpos[$j]-1]"};
		print F @{$samplesizesdigits[$i]}; print F "\r\n"};
	print F "\r\n";
	for (my$i=0; $i<=$#varpos; $i++) {print F '10'};
	close F;
	}
}

print "\n";



