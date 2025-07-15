use Storable;
use EA::GA_Flex_Energy;
use Getopt::Long;
use Data::Dumper;
use strict;

my $usage = "Usage: $0 -cov_model <ranger model file> -energy_model <ranger model file> -init_pop_fn <initial seq population file> -N <population size :default 100> -mut_rate <mutation rate : default 0.1> -survive <survive rate : default 0.05> -iter <iteration number : default 100>\n";

my $SCPT_PATH="../scripts/";		# need to change this path 

my $Seq_length = 17;		# 17 mutation sites  

# create an array of allowed tokens
my @AA = ();
#for ('A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y') {		# excluding Cys 
for ('A',     'D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y') {
	push @AA, $_;
}

# ---- Special Arguments ----
my %args=();
GetOptions("cov_model:s"=>\$args{cov_model},
     			 "energy_model:s"=>\$args{energy_model},
     			 "init_pop_fn:s"=>\$args{init_pop_fn},
     			 "N:i"=>\$args{N},
     			 "mut_rate:i"=>\$args{mut_rate},
     			 "survive:i"=>\$args{survive},
     			 "iter:i"=>\$args{iter},
) or die "Incorrect usage!\n";

my $cov_model = $args{cov_model} or die $usage;
my $energy_model = $args{energy_model} or die $usage;
my $init_pop_fn = $args{init_pop_fn};
my $pop_size = $args{N} || 100;
my $mut_rate = $args{mut_rate} || 0.1;
my $survive_rate = $args{survive} || 0.05;
my $iter = $args{iter} || 100;


my $p = EA::GA_Flex_Energy->new(\&fitness_function,\@AA, $Seq_length);

$p->set_mut_rate($mut_rate);
$p->set_survive_rate($survive_rate);

# initialise the population
my @init_seqs;
if(defined $init_pop_fn){
	open(IN,$init_pop_fn) or die "can't open file:$init_pop_fn\n";
	while(<IN>){
		chomp;
		next if(/^$/);
		push @init_seqs,$_;
	}
	close IN;
	
	$pop_size = @init_seqs; # pop size will be overwritten by the size of seqs 
	$p->init_pop_w_seqary(\@init_seqs);

}elsif(defined $pop_size){
	
	$p->init_pop($pop_size,$Seq_length);

}else{
	
	die $usage;
}

$p->dump_indivs();

printf STDERR "[%5d]\t", 0;	# ??, idx 
$p->dump_best();

my $gen; 

do {
# breed the next generation using crossover and mutation
	$gen = $p->breed();
	printf STDERR "[%5d]\t", $gen-1;	# ??, idx 

# built in data dumper

	$p->dump_best();

# best_fit return the fitness of the best
} while ($p->best_fit() > 0 && $gen <= $iter);

print "\n----------------------------------\n";
print "Cov Cls Model $cov_model\n";
print "Energy Reg Model $energy_model\n";
print "POPULATION $pop_size\n";
print "MutationRate $mut_rate\n";
print "SurviveRate $survive_rate\n";
print "Iterations $iter\n";
print "----------------------------------\n\n";

$p->dump_indivs();

exit(0);


# Now all we need is the fitness function that needs to understand
# the representation of an individual
sub fitness_function {
	
	my ($pop_ref) = @_ ;
	
	my $fn = &write_pop_feat($pop_ref);
	my $score_fn_cov = "$fn.cov"; 
	my $score_fn_energy = "$fn.energy"; 
	
	`Rscript  $SCPT_PATH/predict_regression_model_for_Energy.r $cov_model $fn > $score_fn_cov`;
	`Rscript  $SCPT_PATH/predict_regression_model_for_Energy.r $energy_model $fn > $score_fn_energy`;
	
	my @fit_scores_cov;
	open(IN,$score_fn_cov) or die "can't open file:$_\n";
	while(<IN>){
		chomp;
		my @tmp=split;
		push @fit_scores_cov, $tmp[-1];
	}
	close IN;
	my @fit_scores_energy;
	open(IN,$score_fn_energy) or die "can't open file:$_\n";
	while(<IN>){
		chomp;
		my @tmp=split;
		push @fit_scores_energy, $tmp[-1];
	}
	close IN;

	return (\@fit_scores_cov,\@fit_scores_energy);
}

sub write_pop_feat {
	my ($pop_ref)=@_;
	my $fn1="pop$$.feat1";
	my $fn2="pop$$.feat2";
	open(Out,">$fn1") or die "can't open file:$fn1\n";

	my @sums;
	my $pcnt=1;
	for my $i (@{$pop_ref}){
		my @indiv_tokens = @{$i} ;
		my $s = join "", @indiv_tokens;
		print Out "$s\n";
		$pcnt++;
	}	
	close Out;

	`Rscript $SCPT_PATH/seqfeat_gen.r $fn1 > $fn2`;
	#print STDERR "$fn2 feature file is generated !\n";
	return $fn2;
}

