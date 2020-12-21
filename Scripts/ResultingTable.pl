use strict;
open(CONF, "FitnessShift.Config.txt") or die $!;

my $str1 = <CONF>;
my $str2 = <CONF>;
my $str3 = <CONF>;
my $str4 = <CONF>;
my $str = <CONF>;

my $threshold_line = <CONF>;
chomp($threshold_line);
$threshold_line =~ m/\=/;
my $threshold = "$'";

my $out_line = <CONF>;
chomp($out_line);
$out_line =~ m/\=/;
my $out_file = "$'";
open (OUT1, ">>".$out_file) or die $!;

my %site_aa_subtype_pv1;
my %site_aa_subtype_pv2;
my %site_aa_subtype_min_pv1;
my %site_aa_subtype_min_pv2;

my %site_aan;
my @mafft_sites = <*.site>;
foreach my $ms(@mafft_sites) {
    $ms =~ m/\d+/;
    my $st = $&;
    open (MS, $ms) or die $!;
    my @ms_aas;
    while (<MS>) {
	my $id = $_;
	my $let = <MS>;
	if ($let ~~ @ms_aas) {
	}
	else {
	    push(@ms_aas,$let);
	}
    }
    my $ms_aan = $#ms_aas+1;
    $site_aan{$st} = $ms_aan;
}


#open (SITES, "/mnt/lustre/galkaklink/GAG/GAG.AA_SITES") or die $!;
#my %al_seq;
#my $str1 = <SITES>;
#while (<SITES>) {
#    my $str = $_;
#    my @all = $str =~ m/\w+|\-/g;
#    my $seq = $all[1];
#    my $al = $all[3];
#    $al_seq{$al} = $seq;
#}

my @clades;
my @files = <*p_value.10000>;
foreach my $file (@files) {
    $file =~ m/\d+\.\w+/g;
    my $site_subt = $&;
    $site_subt =~ m/\./;
    my $site = $`;
    my $subt = "$'";

    if ($subt ~~ @clades) {
    }
    else {
	push(@clades,$subt);
    }
    open (IN, $file) or die $!;
    while (<IN>) {
	my $str = $_;
	my @all = $str =~ m/\d+\.\d+|\w+/g;
	my $let = $all[0];
	my $close_pv = $all[1];
	my $min_close_pv = $all[2];
	my $far_pv = $all[3];
	my $min_far_pv = $all[4];
	my $site_aa = $site."_".$let;	
	${$site_aa_subtype_pv1{$site_aa}}{$subt} = $close_pv;
	${$site_aa_subtype_pv2{$site_aa}}{$subt} = $far_pv;
    }
}

print OUT1 "site\tAA_number_in_site\tAA";
foreach my $clade (sort {$a cmp $b} @clades) {
    print OUT1 "\t".$clade."\tAA_z-score_".$clade;
}
print OUT1 "\n";

my @sign_sites;	
foreach my $site_aa (sort {$a <=> $b} keys %site_aa_subtype_pv1) {
    $site_aa =~ m/\_/;
    my $site = $`;
    my $aa = "$'";
    my %subt_pv1 = %{$site_aa_subtype_pv1{$site_aa}};
    my %subt_pv2 = %{$site_aa_subtype_pv2{$site_aa}};
    my @p_values;
    my $flag1 = 0;
    my $flag2 = 0;
    my @close_subt;
    my @far_subt;
    my $aan = $site_aan{$site};
    foreach my $subt (@clades) {
	if (exists $subt_pv1{$subt}) {
	    if ($subt_pv1{$subt} <= $threshold) {
		$flag1 ++;
		push(@close_subt,$subt);
	    }
	    if ($subt_pv2{$subt} <= $threshold) {
		$flag2 ++;
		push(@far_subt,$subt);
	    }
	    push (@p_values,$subt_pv1{$subt});
	    push (@p_values,$subt_pv2{$subt});
	}
	else {
	    push (@p_values,"NA");
	    push (@p_values,"NA");
	}
    }
    if ($flag1 > 0 && $flag2 > 0) {

	print OUT1 $site."\t".$aan."\t".$aa;
       
	foreach my $clade (sort {$a cmp $b} @clades) {
	    my $type = "none";
	    if ($clade ~~ @close_subt) {
		$type = "proximal";
	    }
	    elsif($clade ~~ @far_subt) {
		$type = "distal";
	    }
	    print OUT1 "\t".$type;


	    my $es_file = $site.".".$clade.".z_score.10000";
	    my $z_score = "?";
	    open (ES, $es_file) or die $!;
	    my $names = <ES>;
	    while (<ES>) {
		my $str = $_;
		my @all = $str =~ m/\-\d+\.\d+|\d+\.\d+|\w+/g;
		my $amino_acid = $all[0];
		my $es = $all[$#all];
		if ($amino_acid eq $aa) {
		    $z_score = $es;
		    last;
		}
	    }
	    print OUT1 "\t".$z_score;
	}
	print OUT1 "\n";
    }
}	
