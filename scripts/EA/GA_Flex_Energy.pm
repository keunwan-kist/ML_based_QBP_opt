
# Author: Stefan Kamphausen <mail@skamphausen.de>
# Copyright 2001 Stefan Kamphausen.
# This implements a Simple Somewhat Generalized Genetic Algorithm
# See the bottom of this file for the POD documentation.


####################################################################
##                             LICENSE
####################################################################
# This program is free software; you can redistribute it
# and/or modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
# 02111-1307, USA.

####################################################################
### Original code modified by Keunwan Park (keunwan@kist.re.kr), July 2025
### Modified GA selection strategy, customized scoring function, etc
#####################################################################


package EA::GA_Flex_Energy;

use strict;
use Data::Dumper;

$GA::VERSION='0.1';

sub new {
  my ($class,$fitness,$tokenref,$seq_len) = @_;
  my $self = {};
  # An array which contains all the individuals as arrays of tokens
  $self->{pop} = ();
  # A Hash which contains all the fitness-values adressed by
  # the concatenated tokens
  $self->{fitvals_cov} = ();
  $self->{fitvals_energy} = ();
  $self->{fitvals_combined} = ();
  # the 'alphabet' of allowed symbols
  $self->{tokens} = $tokenref;
  # user provided fitness function, gets array of tokens as arg
  $self->{fitness} = $fitness;
  # default value for mutation probability; may be overridden
  $self->{mut_prob} = 0.1;
  # a counter for the generation
  $self->{generation} = 0;
  # Survived population percentage 
  $self->{survive_rate} = 0.05;
  # Pop entry seq length  
  $self->{seq_length} = $seq_len;
  bless $self, $class;
  return $self;
}

sub init_pop {
  my ($self,$size,$seq_len) = @_;
  for (my $i=0;$i<$size;$i++) {	# population size
	for (my $j=0;$j<$seq_len;$j++) {	# sequence length
	  my $rtok = $self->random_token();
	  #print "RTOK: $rtok\n";
	  push @{$self->{pop}[$i]}, $rtok;
	}
	#my $s = join("",@{$self->{pop}[$i]});
	#print "$s\n";
  }
  $self->calculate_fitness();
  $self->sort_pop();
  ++$self->{generation};
}
sub init_pop_w_seqary {
  my ($self,$seq_ary) = @_;
  my $size = keys @{$seq_ary};
  my $seq_len = length $seq_ary->[0];
  for(my $i=0;$i<$size;$i++) {	# population size
	my $s = $seq_ary->[$i];
	my @ss=split(//,$s);
	@{$self->{pop}[$i]} = @ss;
	#my $s = join("",@{$self->{pop}[$i]});
	#print "$s\n";
  }
  $self->calculate_fitness();
  $self->sort_pop();
  ++$self->{generation};
}

sub show_pop {
  my ($self) = @_;
  my @tmp_pop = @{$self->{pop}};

  for my $p (@tmp_pop){
	my @aa = @{$p};
	my $aa_str = join("",@aa);
	print "$aa_str\n";
  }
}

sub sort_pop {
  my ($self) = @_;
  @{$self->{pop}} = sort {$self->{fitvals_combined}{join "",@{$b}} <=> $self->{fitvals_combined}{join "",@{$a}}} @{$self->{pop}};
}

sub calculate_fitness {
  my ($self) = @_;
  
  my $gen = $self->{generation};
  %{$self->{fitvals_cov}} = ();
  %{$self->{fitvals_energy}} = ();
  %{$self->{fitvals_combined}} = ();
  my ($fits_ref_cov,$fits_ref_energy) = &{$self->{fitness}}(\@{$self->{pop}});	# treat all elements, renew $self->{fitvals}{seq}
 
  for my $i (@{$self->{pop}}) {
	my $s = join "",@{$i};
	my $f1=shift @{$fits_ref_cov};   
	$self->{fitvals_cov}{$s} = $f1;
	my $f2=shift @{$fits_ref_energy};   
	$self->{fitvals_energy}{$s} = $f2;

	# my $combined_f = -1*($f2/850);  # 100:0, -850 WT ddg, higher better
	# my $combined_f = -1*($f1/26)*0.8 + -1*($f2/850)*0.2 ; # 80:20
	# my $combined_f = -1*($f1/26)*0.5 + -1*($f2/850)*0.5 ; # 50:50
	# my $combined_f = -1*($f1/26)*0.2 + -1*($f2/850)*0.8 ; # 20:80
	
	my $combined_f = -1*($f1/26);	# 0:100, -26 WT ddg, higher better  
	
	$self->{fitvals_combined}{$s} = $combined_f;	
	print STDERR "$gen Indiv:$s $combined_f $f1 $f2\n";
  }
}

### Return some values
sub best_fit {
  my $self = shift;
  return $self->{fitvals_combined}{join("",@{$self->{pop}[0]})};
}
sub generation {
  my $self = shift;
  return $self->{generation};
}
sub random_token {
  my ($self) = @_;
  my $max = scalar(@{$self->{tokens}});
  
  my $ran = int rand $max;
  my $tok = @{$self->{tokens}}[$ran];
  #print "Ran: $ran, Max: $max ";
  #print "Token: $tok\n";
  return $tok;
}

# combines mutation and crossover
sub breed {
  my $self = shift;
  my $opt_mutation_rate = shift;

  my ($p1,$p2,$c,$i);
  # mutation
  my $p_mut = $opt_mutation_rate || $self->{mut_prob};
  # prepare for roulette wheel
  my $a_sum = 0;
  my @fit = ();
  foreach (@{$self->{pop}}) {
	push @fit, $self->{fitvals_combined}{join "",@{$_}};
  }

  @fit = sort {$b <=> $a} @fit;

  #print "BREED: Fitness\n",Dumper(@fit),"\n\n";
  my $length = scalar(@{$self->{pop}});
  
  my @new_pop = ();
  my @elit_pop = ();	
  
  # N % survive
  my $survive_n = int( $length * $self->{survive_rate} ); # at least one 
  
  for(my $i=0;$i<=$survive_n;$i++){
	push @elit_pop, @{$self->{pop}}[$i];
	push @new_pop, @{$self->{pop}}[$i];
  }
  my @elit_fit=@fit[0..$survive_n];
  #print STDERR "$survive_n survived : @elit_pop : @elit_fit \n";
  
  # mutations to all population
  $self->mutate($p_mut,$survive_n);	# whole mutations 
 
  # Choose Parents
  for ($i=$survive_n+1;$i<$length;$i++) {
	$p1 = rwheel(\@fit);
	$p2 = rwheel(\@fit);
	#print "P1: $p1 P2: $p2\n";
	
	my $new_gene;
	my $iter=0;
	while($iter<100){
		$new_gene = $self->crossover_mut($p1,$p2,$p_mut);
		my $s = join("",@{$new_gene});
		if(defined $self->{fitvals_combined}{$s}){
			# generate redundant gene, try again for N-times
			#print STDERR "$s is redundant... try again..\n";
		}else{
			# uniq gene
			last;
		}
		$iter++;
	}
	push @new_pop, $new_gene;
  }

  @{$self->{pop}} = @new_pop;

  $self->calculate_fitness();
  $self->sort_pop();
  return ++$self->{generation};
}

sub mutate {
  my ($self,$rate,$survive_n) = @_;
  my ($ran,$i,$t);
  
  my $the_rate = $rate || $self->{mut_prob};
  my $cnt=0;	 
  foreach $i (@{$self->{pop}}) {
	$cnt++;
	next if($cnt <= $survive_n+1 );	# starting from 1
	for ($t=0;$t<scalar(@{$i});$t++) {
	  $ran = rand();
	  if ($ran < $the_rate) {
		#print "MUTATE $cnt-th seq $t-th pos!\n";
		@{$i}[$t] = $self->random_token();
	  }
	}
  }
}

sub crossover_mut {
  my $self = shift;
  my $p1 = shift;
  my $p2 = shift;
  my $opt_mutation_rate = shift;
  
  my ($ran,$t,$new_size);
  $ran = rand();
  my $pmut = $opt_mutation_rate || $self->{mut_prob};
  
  # 50:50 for the size of the new one 
  my @new = @{$self->{pop}[$p1]};
  my @new1 = @{$self->{pop}[$p1]};
  my @new2 = @{$self->{pop}[$p2]};
  my $new_size = @new;	
  
  # mutation	
  for($t=0;$t<@new;$t++) {
	$ran = rand();
	if ($ran < $pmut) {
	  $new1[$t] = $self->random_token();
	}
	$ran = rand();
	if ($ran < $pmut) {
	  $new2[$t] = $self->random_token();
	}
  }
  
  my $pp1 = 0.5;
  for($t=0;$t<$new_size;$t++) {
	$ran = rand();
	# 50:50 to take gene from p1 or p2 unless mutation
	if($ran < $pp1) {
	  $new[$t] = $new1[$t];	  
	} else {
	  $new[$t] = $new2[$t];	  
	}
  }
  return \@new;
}

sub crossover {
  my $self = shift;
  my $p1 = shift;
  my $p2 = shift;
  
  my ($ran,$t,$new_size);
  $ran = rand();
  # 50:50 for the size of the new one 
  my @new = ();
  if ($ran < 0.5) {
	$new_size = scalar(@{$self->{pop}[$p1]})
  } else {
	$new_size = scalar(@{$self->{pop}[$p2]})	
  }
  for($t=0;$t<$new_size;$t++) {
	
	$ran = rand();
	# 50:50 to take gene from p1 or p2 unless mutation
	if ($ran < 0.5) {
	  $new[$t] = @{@{$self->{pop}}[$p1]}[$t];	  
	} else {
	  $new[$t] = @{@{$self->{pop}}[$p2]}[$t];	  
	}
  }
  return \@new;
}


### Print-Outs
sub dump_indivs {
  my $self = shift;
  my $i;
  my $len = scalar(@{$self->{pop}});
  for ($i=0;$i<$len;$i++) {
	my $s = join("",@{$self->{pop}[$i]});
	printf "%4d %s %5f %5f %5f\n", $i, $s, $self->{fitvals_combined}{$s}, $self->{fitvals_cov}{$s}, $self->{fitvals_energy}{$s};
  }
}

sub calc_avg_fit {
  my $self = shift;
  my $i;
  my $len = scalar(@{$self->{pop}});
  my $avg=0;
  for ($i=0;$i<$len;$i++) {
	my $s = join("",@{$self->{pop}[$i]});
	$avg += $self->{fitvals_combined}{$s};
  }
  $avg = sprintf("%3.4f",$avg/$len);
  return $avg;
}


sub dump_best {
  my $self = shift;
  my $s = join("",@{$self->{pop}[0]});
  my $avg = $self->calc_avg_fit();
  printf STDERR "%s %5f %5f %5f Avg %5f\n", $s, $self->{fitvals_combined}{$s}, $self->{fitvals_cov}{$s}, $self->{fitvals_energy}{$s},$avg;
}

sub set_mut_rate {
  my ($self,$mut_rate) = @_;
  $self->{mut_prob} = $mut_rate;
}
sub set_survive_rate {
  my ($self,$survive_rate) = @_;
  $self->{survive_rate} = $survive_rate;
}
sub set_seq_len {
  my ($self,$slen) = @_;
  $self->{seq_length} = $slen;
}

### Random
sub rwheel {
  # random element of an array according to it's value
  # aka roulette wheel
  my ($a_ref) = @_;
  my @arr = @{$a_ref};
  
  my $a_sum = 0;
  foreach (@arr){
	$a_sum+= $_
  }
  
  my $sum = 0;
  #print "RWHEEL: length = ",scalar(@arr),"\n";
  #print "RWHEEL ARRAY: ",join(" ",@arr),"\n";
  my $ran = rand $a_sum;
  #print "RWHEEL: RAN $ran < $a_sum\n";
  for (my $i=0;$i<scalar(@arr);$i++) {
	$sum += $arr[$i];
	#print "\tSUM: $sum \$arr[$i] = $arr[$i]\n";
	if ($sum > $ran ) {
	  return $i;
	}
  }
  die "ARGH! I never should have reached this point!\n";
}

1;

__END__
############################################################
#                           DOCS                           #
############################################################

=head1 NAME

EA::GA - a general genetic algorithm library

=head1 SYNOPSIS

    # This is a little example

    use EA::GA;

    # evolve a string that matches this target
    $target = "Hello_World";
    $len = length $target;

    # create an array of allowed tokens
    @token = ();
    for ('a'..'z') {
      push @token, $_;
    }
    for ('A'..'Z') {
      push @token, $_;
    }
    push @token, "_";

    # New GA object that sets the alphabet and the
    # fitness function
    $p = EA::GA->new(\&fitness_function,\@token);

    # initialise the population
    $p->init_pop(100,$len);

    do {
      # breed the next generation using crossover and mutation
      $gen = $p->breed();
      printf "[%5d] ", $gen;
      # built in data dumper
      $p->dump_best();
      # best_fit return the fitness of the best
    } while ($p->best_fit() > 0 && $gen < 2000);

    $p->dump_best();
    exit(0);

    # Now all we need is the fitness function that needs to understand
    # the representation of an individual

    sub fitness_function {
      my @indiv_tokens = @_ ;
      # Representation
      my $s1 = join "", @indiv_tokens;
      my $sum = 0;
      my $f;
      for($f=0;$f<$len;$f++) {
        my $z1=substr($s1,$f,1);
        my $z2=substr($target,$f,1);
        my $a=(ord($z1)-ord($z2))*(ord($z1)-ord($z2));
        $sum +=$a;
      }
      return $sum;
    }

=pod

=head1 DESCRIPTION

C<EA::GA> implements a (hopefully) generalized genetic algorithm.
It does this by using an array of allowed tokens as individuals.
The user has to provide a fitness function. There the actual
representation is implemented. If you got a string of chars it is
quite easy: simply join them. If you want to have real numbers you
should probably use a bitwise representation and calculate the
real values in your fitness function.

=head2 The Easy Way

The easy setup is pretty easy. With

   $p = EA::GA->new(\&fitness_function,\@token);

you create a new GA object which knows all the allowed tokens and
how to calculate the fitness of an individual.
Then use

   $p->init_pop($pop_size,$length_of_individual);

to initialise a random population of I<$pop_size> individuals, each
of length I<$length_of_individual>. I do not know how to make them
of variable length right now.

The main thing to do now is use the simplified C<breed()>-method

   $gen = $p->breed();

You can give an optional argument to the C<breed> method which will
be interpreted as the mutation probabiliy. This method combines
mutation and crossover (for each token there is a decision from
which parent to take the token) and returns the number of the
generation.

=head2 The Detailed Way

There are methods that provide mutation, crossover and other
functionality and can be called directly in case you do not want to
use the built in C<breed()> method. These and other methods will soon
be listed in alphabetical order. Right before that again the note that
you probably don't need this.

=over 4

=item best_fit()

Returns the fitness of the best individual of the whole population
if the population is sorted (actually returns the first element of
the internal population array).

=item calculate_fitness()

Updates the (internal) fitness values by calling the user provided
fitness function for each individual.

=item crossover($p1,$p2)

Does a simple crossover schema. All individuals are internally
represented as an array of tokens. This crossover needs the numbers of
two parents (I<$p1> and I<$p2>), usually drawn using the Roulette
Wheel technique. For each token of the offspring there is a
fifty:fifty decision whether to take from parent one or parent two.

=item crossover_mut($p1,$p2,$optional_mutation_prob)

Almost the same as C<crossover()> just that there is a little
probabiliy that a new random token is used instead of on of the
parents.

=item dump_best()

This prints the best individual to stdout in a somewhat reasonable
way.

=item dump_indivs()

Prints the whole population including their fitness values.

=item generation()

Returns the number of the current generation.

=item mutate($optional_mutation_prob)

Performs a mutation on the whole generation.

=item sort_pop()

Whenever a new population has been created and the fitness values have
been calculated it is necessary to sort the population. Some routines
rely on that.

=item random_token()

Return a random token from the user provided alphabet of allowed
tokens. 

=back

=head1 AUTHORS

Stefan Kamphausen I<E<lt>mail@skamphausen.deE<gt>>
I<http://www.skamphausen.de/software>

=cut
