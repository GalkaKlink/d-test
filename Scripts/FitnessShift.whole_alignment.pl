#This script does the same thing that FitnessShift.pl does, but this one is suitable in case that there is no possibilities for parallel computation for protein sites, and you want to analyse each site iteratively. 
use strict;
use Bio::TreeIO;

open(CONF, "FitnessShift.Config.txt") or die $!;
my %subtype_id;

my $file = <CONF>;
chomp($file);

my $tree_line = <CONF>;
chomp($tree_line);
$tree_line =~ m/\=/;
my $tree_file = "$'";

my $root_line = <CONF>;
chomp($root_line);
$root_line =~ m/\=/;
my $rt = "$'";

my $clade_line = <CONF>;
$clade_line =~ m/\=/;
my $clade_line1 = "$'";
my @clades = $clade_line1 =~ m/\w+/g;

my $id_line = <CONF>;
$id_line =~ m/\=/;
my $id_line1 = "$'";
my @ids = $id_line1 =~ m/\w+/g;

if ($#clades != $#ids) {
    die "Number of focal_clades must be equal to number of focal_nodes!";
}

foreach my $k(0..$#clades) {
    my $clade = $clades[$k];
    my $id = $ids[$k];
    $subtype_id{$clade} = $id;
    print $clade."\t".$id."\n";
}


my $in = Bio::TreeIO -> new(-file => $tree_file,
			    -format => 'newick');
my $tree = $in -> next_tree;

open(IN,$file) or die $!;

my $id = <IN>;
my $seq = <IN>;
chomp($seq);
my $lenseq = length($seq);

seek (IN,0,0);
    
my %hash = ();
my $k;
my $v;
while (<IN>)  {
    
    $k = $_;
    chomp($k);
    $k =~ s/\>//;

    $v = <IN>;
    chomp($v);
    $hash{$k} = $v;
}

my @leafnodes = $tree->get_leaf_nodes;
my @leafs;
my @neschet;

#########################Next step is needed to exclude terminal branches that are too long (>10 median values of branch lengths) from the analysis. If you do not need this, you can skip code until next "###"
my %leaf_length;
my @brlength;
foreach my $leafnode(@leafnodes)  {
    my $leaf = $leafnode->id;
    push(@leafs,$leaf);
    my $bl = $leafnode -> branch_length;
    $leaf_length{$leaf}=$bl;
   
    push(@brlength,$bl);
}

my @sorted = sort {$a<=>$b} @brlength;

my $middle = ($#sorted+1)/2;
$middle =~ m/\d+/;
my $mid = $&;
my $median = $sorted[$mid-1];
my $tenmed = $median*50;

foreach my $leaf(keys %leaf_length)  {
    my $bl = $leaf_length{$leaf};
    if ($bl >= $tenmed) {
	push(@neschet,$leaf);
    }    	
}		

my $tbl = $tree->total_branch_length;
foreach my $ns_id(@neschet)  {
    my $ns = $tree->find_node(-id=>$ns_id);
    my $ns_bl = $ns->branch_length;
    $tbl = $tbl-$ns_bl;
 }
##################################

my %level;  

my $m = 0;
my $q = 0;
push (@{$level{$m}},$rt);
my $r_id = $rt;

lev ($r_id,$q); #This function makes groups of nodes which belong to the same level from the route (route has level 0, its decsendants have level 1, etc)


my $n = 1;
foreach my $n(1..$lenseq) { #this is needed if parallel calculationes are not possible and you want to give the whole alignment to a program instead of each column in a separate file
    my $site = $n;
    my %from_to_nodes;
    my %fromB;
    my %node_aa; 
    
    my $root = $rt;
    my $root_seq = $hash{$root};
    my $root_aa = substr($root_seq,$n-1,1);
    
    until ($m == -1)  {
	
    my @nds = @{$level{$m}};
    
    foreach my $n_id(@nds)  {
	
	if ($n_id ~~ @leafs)  {
	    
	}
	
	else  {
	    my $node = $tree->find_node(-id=>$n_id);
	    
	    my @des = $node->each_Descendent;
	    
	    foreach my $des(@des) {
		my $n1 = $des->id;
		#    print "des\t".$n1."\n";
		
		if (($n1 ~~ @neschet) || ($n1 == undef)) {
		}
		
		else {
		    my $bl = $des->branch_length;
		    #	print "bl\t".$bl."\n";			
			
		    my $ancc = $des->ancestor;
		    my $n2 = $ancc -> id;
		    
		    my $s1 = $hash{$n1};
		    my $s2 = $hash{$n2};
		    
		    my $aa1 = substr($s1,$n-1,1);
		    my $aa2 = substr($s2,$n-1,1);
		    #	print "aas\t".$aa1."\t".$aa2."\n";
		    
		    if ($aa1 =~ /\w/ && $aa2 =~ /\w/ && $aa1 ne '-' && $aa2 ne '-')  {
			
			if ($aa1 ne $aa2 && $aa1 ne 'X' && $aa2 ne 'X')  {
			    
			    push(@{${$from_to_nodes{$aa2}}{$aa1}},$n1);
			    push(@{$fromB{$aa2}},$n1);
			    $node_aa{$n1} = $aa1; 
			    
			}
		    }
		    
		}
	    }	
	    
	}	    
    }
    $m=$m-1;
    }
    
    foreach my $from (keys %from_to_nodes)  {
	print "from\t".$from."\n";
	my %to_nodes = %{$from_to_nodes{$from}};
	my @fB = @{$fromB{$from}};
	
	if (keys %to_nodes > 0)  {
	    my %to_dists;
	    
	    foreach my $n2 (@fB)  {
		foreach my $subtype(keys %subtype_id) {
		    my $hs_id_A = $subtype_id{$subtype};
		    my $n11 = $hs_id_A;
		    my $krugl_distA = 0.000;
		    
		    if ($n11 != $n2) {
			my $node1 = $tree->find_node(-id=>$n11);
			my $node2 = $tree->find_node(-id=>$n2);
			my @nodes = ($node1,$node2);
			my $lca = $tree->get_lca(-nodes=>\@nodes);
			my $lca_id = $lca->id;
			if($lca_id != $n2)  {
			    my $distance = $tree->distance(-nodes=>[$node1,$node2]);
			    my $bl1 = $node1->branch_length;
			    my $bl2 = $node2->branch_length;
			    my $bls = ($bl1+$bl2)/2;
			    $distance = $distance-$bls;
			    
			    if ($distance =~ m/e/|$distance =~ m/0\.000/ && $distance<1)  {
				$krugl_distA = 0.001;
			    }
			    elsif ($distance =~ m/\./)  {
				$distance =~ m/\./;
				my $int = $`;
				my $dec = "$'";
				my @numb = $dec =~ m/\d/g;
				my $one = $numb[0];
				my $two = $numb[1];
				my $three = $numb[2];
				my $four = $numb[3];
				if ($four>=5) {
				    $three+=1;
				    $four = 0;
				    if ($three>=5) {
					$two+=1;
					$three = 0;
					if ($two>=5)  {
					    $one +=1;
					    $two = 0;
					    if ($one == 10) {
						$int+=1;
						$one = 0;
					    }
					}
				    }
				}
				else  {
				}
				
				$krugl_distA = $int.".".$one.$two.$three;
			    }
			    else {
				$krugl_distA = $distance.".000";
			    }
			}
		    }
		    
		    my $aa = $node_aa{$n2};
		    push(@{${$to_dists{$subtype}}{$aa}},$krugl_distA);
		}
	    }
	    
	    foreach my $subtype(keys %to_dists) {
		my $id = $subtype_id{$subtype};
		
		open (OUT1, ">>$site.dists_from".$subtype."_".$id) or die $!;
		my %to_distsA = %{$to_dists{$subtype}};
		
		if (keys %to_distsA > 0) {
		    foreach my $aa(keys %to_distsA) {
			print OUT1 $site."\t".$from."\t".$aa;
			my @distsA = @{$to_distsA{$aa}};
			foreach my $dist(@distsA) {
			    print OUT1 "\t".$dist;
			}	
			print OUT1 "\n";
		    }
		}
		
	    }
	}
    }

###################NEXT STEP: calculation of p-value, z-scores and prevalences
    my @files = <$site.dists_from*>;
    
    foreach my $file(@files) {
	$file =~ m/dists_from/;
	my $clade_id = "$'";
	$clade_id =~ m/\_/;
	my $clade = $`;
	my $id = "$'";
	
	open(OUT2, ">>".$site.".".$clade_id.".p_value.10000") or die $!;
	open(OUT3, ">>".$site.".".$clade_id.".z_score.10000") or die $!;
	print OUT3 "AA\tmean_dist\tH0_mean_dist\tH0_std_dev\tz-score\n";
	open (IN, $file) or die $!;
	my %from_dists;
	my %to_dists;
	my %to_from_numb;
	while (<IN>) {
	    my $str = $_;
	    my @all = $str =~ m/\w+/g;
	my $from = $all[1];
	    my $to = $all[2];
	    my @dists = $str =~ m/\d+\.\d/g;
	    my $chisl;
	    my $znam;
	    foreach my $dist(@dists) {
		${$to_from_numb{$to}}{$from} += 1;
	    push(@{$from_dists{$from}},$dist);
		push(@{$to_dists{$to}},$dist);
	    }
	}
	
	foreach my $to (keys %to_dists) {
	    my @dists = @{$to_dists{$to}};
	    my $to_number = $#dists+1;
	    #print $to."\t".$to_number."\n";
	    
	    my %from_numb = %{$to_from_numb{$to}};
	    my $from_number = 0;
	    foreach my $from(keys %from_numb) {
		my @fd = @{$from_dists{$from}};
		my $numb = $#fd+1;
		$from_number += $numb;
		#print $from."\t".$numb."\n";
	    }
	    
	    
	    my $mean;
	    foreach my $dist(@dists) {
		$mean += $dist;
	    }
	    my $sum_dist = $mean;
	    $mean=$mean/($#dists+1);
	    
	    my @sim_sum_dists;
	    #my @sim_vars;
	    until ($#sim_sum_dists == 9999) {
		my $sim_sum_dist;
		my @sim_dists;
		my %from_numb = %{$to_from_numb{$to}}; 
		foreach my $from (keys %from_numb) {
		    my $numb = $from_numb{$from};
		    my @dists= @{$from_dists{$from}};
		    my $chisl;
		    my @used;
		    my $dists_number = $#dists+1;
		    until($#used == $numb-1) {
			my $rand = int(rand($dists_number));
			if ($rand~~ @used) {
			}
			else {
			    $sim_sum_dist+=$dists[$rand];
			    push(@used,$rand);
			    push(@sim_dists,$dists[$rand]);
			}
		    }
		}
		push(@sim_sum_dists,$sim_sum_dist);
		my $sim_mean = $sim_sum_dist/($#sim_dists+1);   
    }
	    
	    my @sort_sim_sum_dists = sort {$a <=> $b} @sim_sum_dists;
	    
	    
	    #for p-values:
	    my $min_dist = $sort_sim_sum_dists[0];
	    my $min_pv = 1;
	    foreach my $n(0..$#sort_sim_sum_dists) {
		if ($sort_sim_sum_dists[$n]>$min_dist) {
		    $min_pv = $n/10000;
		    last;
		}
	    }
	    my $pv = 1;
	    foreach my $n(0..$#sort_sim_sum_dists) {
		if ($sort_sim_sum_dists[$n]>$sum_dist) {
		    $pv = $n/10000;
		    last;
		}
	    }
	    
	    
	    
	    my @sort_sim_sum_dists_rev = sort {$b <=> $a} @sim_sum_dists;
	    my $max_dist = $sort_sim_sum_dists_rev[0];
	    my $min_pv_rev = 1;
	    foreach my $n(0..$#sort_sim_sum_dists_rev) {
		if ($sort_sim_sum_dists_rev[$n]<$max_dist) {
		    $min_pv_rev = $n/10000;
		    last;
		}
	    }
	    my $pv_rev = 1;
	    foreach my $n(0..$#sort_sim_sum_dists_rev) {
		if ($sort_sim_sum_dists_rev[$n]<$sum_dist) {
		    $pv_rev = $n/10000;
		    last;
		}
	    }
	    
	    if ($pv == 0) {
		$pv = $min_pv;
	    }
	    if ($pv_rev == 0) {
		$pv_rev = $min_pv_rev;
	    }
	    
	    print OUT2 $to."\t".$pv."\t".$min_pv."\t".$pv_rev."\t".$min_pv_rev."\n";
	    
	    #for z-scores:
	    my $chisl = 0;
	    my $znam = 0;
	    foreach my $dist(@sort_sim_sum_dists) {
		$chisl+=$dist;
		$znam+=1;
	    }
	    my $distr_mean_sum_dist = $chisl/$znam;
	    my $eff_chisl = $sum_dist/$to_number;
	    my $eff_znam =$distr_mean_sum_dist/$to_number;
	    my $eff_size = $eff_chisl/$eff_znam;
	    
	    my $mean_H0 = &average(\@sort_sim_sum_dists);
	    my $std_H0 = &stdev(\@sort_sim_sum_dists);
	    my $z_score = "NA";
	    if ($std_H0 > 0) {
		$z_score = ($sum_dist-$mean_H0)/$std_H0;
	    }
	    print OUT3 $to."\t".$sum_dist/$to_number."\t".$mean_H0/$to_number."\t".$std_H0/$to_number."\t".$z_score."\n";
	}
    }
}




sub lev  {
    
    my ($par_id,$n) = @_;
        
    $n++;
#    print $par_id."\n";
    my $par = $tree->find_node(-id => $par_id);

  CHILD: for my $child ($par -> each_Descendent)  {
      
      my $ch_id = $child->id;
      push (@{$level{$n}},$ch_id);

      if ($n > $m)  {
	  
	  $m = $n;
      }

      if($ch_id ~~ @leafs) {	  
      }
      
      else {
	  lev ($ch_id,$n);	  
      }
      next CHILD;
  }
    
}


sub average{
    my($data) = @_;
    if (not @$data) {
        die("Empty arrayn");
    }
    my $total = 0;
    foreach (@$data) {
        $total += $_;
    }
    my $average = $total / @$data;
    return $average;
}


sub stdev{
    my($data) = @_;
    if(@$data == 1){
        return 0;
    }
    my $average = &average($data);
    my $sqtotal = 0;
    foreach(@$data) {
        $sqtotal += ($average-$_) ** 2;
    }
    my $std = ($sqtotal / (@$data-1)) ** 0.5;
    return $std;
}

