=head1 NAME

Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranscriptUtils - utilities for transcript objects

=head1 SYNOPSIS

  use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranscriptUtils qw(clone_Transcript);

  or 

  use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranscriptUtils 

  to get all methods

=head1 DESCRIPTION

All methods in this class should take a Bio::EnsEMBL::Transcript
object as their first argument.

The methods provided should carry out some standard 
functionality for said objects such as printing info, and 
cloning and checking phase consistency or splice sites etc

=head1 CONTACT

please send any questions to ensembl-dev@ebi.ac.uk

=head1 METHODS

the rest of the documention details the exported static
class methods

=cut

package Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranscriptUtils;

use strict;
use warnings;
use Exporter;

use Bio::EnsEMBL::Utils::Exception qw(verbose throw warning
                                      stack_trace_dump);
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranslationUtils qw(print_Translation clone_Translation);
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::ExonUtils qw(print_Exon clone_Exon Exon_info 
                                                                exon_length_less_than_maximum);
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::IntronUtils qw(intron_length_less_than_maximum get_splice_sites);
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils qw(coord_string id);
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::EvidenceUtils qw (print_Evidence clone_Evidence);
use Bio::EnsEMBL::Transcript;
use Bio::EnsEMBL::Translation;
use Bio::EnsEMBL::Analysis::Tools::Logger qw(logger_verbosity
                                             logger_info
                                             logger_warning);
use Bio::EnsEMBL::Analysis;
use Bio::EnsEMBL::Analysis::Runnable::ProteinAnnotation::Seg;

use vars qw (@ISA @EXPORT);

@ISA = qw(Exporter);
@EXPORT = qw(
             print_Transcript
             clone_Transcript
             print_Transcript_evidence
             Transcript_info
             are_strands_consistent
             lies_inside_of_slice
             exon_lengths_all_less_than_maximum
             intron_lengths_all_less_than_maximum
             are_phases_consistent
             is_not_folded
             low_complexity_less_than_maximum
             has_no_unwanted_evidence
             is_spliced 
             are_splice_sites_canonical
             count_non_canonical_splice_sites
             exonic_proportion
             coding_coverage
             list_evidence
             split_Transcript
             replace_stops_with_introns
            );




=head2 print_Transcript

  Arg [1]   : Bio::EnsEMBL::Transcript
  Arg [2]   : string, this should be a string or spaces or tabs
  to indent the printed string
  Function  : print information about the transcript and its
  children objects, using indent to make the format readable
  Returntype: n/a
  Exceptions: none
  Example   : print_Transcript($transcript);

=cut


sub print_Transcript{
  my ($transcript, $indent) = @_;
  $indent = '' if(!$indent);
  print Transcript_info($transcript, $indent)."\n";
  my $translation_indent = $indent."\t";
  print_Translation($transcript, $translation_indent);
  foreach my $exon(@{$transcript->get_all_Exons}){
    my $exon_indent = $translation_indent."\t";
    print_Exon($exon, $exon_indent);
  }
}


=head2 Transcript_info

  Arg [1]   : Bio::EnsEMBL::Transcript
  Arg [2]   : string, indent
  Function  : return string of info about the transcript
  Returntype: 
  Exceptions: none
  Example   : print_just_Transcript($transcript);

=cut



sub Transcript_info{
  my ($transcript, $indent) = @_;
  $indent = '' if(!$indent);
  my $coord_string = coord_string($transcript);
  my $id = id($transcript);
  return $indent."TRANSCRIPT: ".$id." ".$coord_string."\n";
}

=head2 print_Transcript_evidence

  Arg [1]   : Bio::EnsEMBL::Transcript
  Arg [2]   : string, an indent
  Function  : print the transcripts supporting evidence
  Returntype: n/a
  Exceptions: none
  Example   : print_Transcript_evidence($transcript);

=cut



sub print_Transcript_evidence{
  my ($transcript, $indent) = @_;
  print $indent."TRANSCRIPT EVIDENCE:\n";
  foreach my $evidence(@{$transcript->get_all_supporting_features}){
    my $evidence_indent = $indent."\t";
    print_Evidence($evidence, $evidence_indent);
  }
}



=head2 clone_Transcript

  Arg [1]   : Bio::EnsEMBL::Transcript
  Function  : produce a new copy of the transcript object passed
  in so it can be altered without impact on the original objects
  the only bit it doesnt keep is the adaptor so the cloned 
  object can be stored
  Returntype: Bio::EnsEMBL::Transcript
  Exceptions: none 
  Example   : my $newtranscript = clone_Transcript($transcript);

=cut



sub clone_Transcript{
  my ($transcript) = @_;
  my $newtranscript = new Bio::EnsEMBL::Transcript();
  foreach my $exon(@{$transcript->get_all_Exons}){
    my $newexon = clone_Exon($exon);
    $newtranscript->add_Exon($newexon);
  }
  foreach my $sf(@{$transcript->get_all_supporting_features}){
    my $newsf = clone_Evidence($sf);
    $newtranscript->add_supporting_features($newsf);
  }
  my $newtranslation = clone_Translation($transcript, 
                                         $newtranscript);
  $newtranscript->translation($newtranslation);
  my $attribs = $transcript->get_all_Attributes();
  $newtranscript->add_Attributes(@$attribs);
  $newtranscript->slice($transcript->slice);
  $newtranscript->dbID($transcript->dbID);
  return $newtranscript;
}



=head2 are_strands_consistent

  Arg [1]   : Bio::EnsEMBL::Transcript
  Function  : checks if strand is consistent between 
  transcript and first exon and in multiexon genes between 
  all exons
  Returntype: boolean, 1 if pass 0 if fail (ie strands 
                                                inconsistent)
  Exceptions: none
  Example   : throw("Strands not consistent") 
  if(!are_strands_consistent($transcript));

=cut



sub are_strands_consistent{
  my ($transcript) = @_;
  my $exons = $transcript->get_all_Exons;
  if($exons->[0]->strand != $transcript->strand){
    logger_warning("Strands are inconsistent between the ".
                   "first exon and the transcript for ".
                   id($transcript));
    return 0;
  }
  if(@$exons >= 2){
    for(my $i = 1;$i < @$exons;$i++){
      if($exons->[$i]->strand != $exons->[$i-1]->strand){
        logger_warning("Strands are inconsistent between ".
                       "exon $i exon and exon ".($i-1)." for ".
                       id($transcript));
        return 0;
      }
    }
  }
  return 1;
}



=head2 lies_inside_of_slice

  Arg [1]   : Bio::EnsEMBL::Transcript
  Arg [2]   : Bio::EnsEMBL::Slice
  Function  : ensures the transcript within the slice, 
  completely on the lower end, it can overhang the upper end
  Returntype: boolean, 1 for pass 0 for fail ie(lies outside
                                                    of slice or
                                                    across lower 
                                                    boundary)
  Exceptions: none
  Example   : 

=cut


sub lies_inside_of_slice{
  my ($transcript, $slice) = @_;
  if($transcript->start > $slice->length || 
     $transcript->end < 1){
    logger_warning(id($transcript)." lies off edge if slice ".
                   $slice->name);
    return 0;
  }
  if($transcript->start < 1 && $transcript->end > 1){
    logger_warning(id($transcript)." lies over lower boundary".
                   " of slice ".$slice->name);
    return 0;
  }
  return 1;
}



=head2 exon_lengths_all_less_than_maximum

  Arg [1]   : Bio::EnsEMBL::Transcript
  Arg [2]   : int, max length
  Function  : checks if any of the exons of given transcript
  are longer than specified max length
  Returntype: boolean, 1 for pass, 0 for fail (ie exon beyond
                                                   max length)
  Exceptions: none
  Example   : 

=cut



sub exon_lengths_all_less_than_maximum{
  my ($transcript, $max_length) = @_;
  foreach my $exon(@{$transcript->get_all_Exons}){
    if(!exon_length_less_than_maximum($exon, $max_length)){
      logger_info("Transcript ".id($transcript)." has ".
                  "exon longer than ".$max_length);
      return 0;
    }
  }
  return 1;
}

=head2 intron_lengths_all_less_than_maximum

  Arg [1]   : Bio::EnsEMBL::Transcript
  Arg [2]   : int, max length
  Function  : checks if any of the introns of given transcript
  are longer than specified max length
  Returntype: boolean, 1 for pass, 0 for fail 
  (ie intron beyond max length)
  Exceptions: none
  Example   : 

=cut



sub intron_lengths_all_less_than_maximum{
  my ($transcript, $max_length) = @_;
  foreach my $intron(@{$transcript->get_all_Introns}){
    if(!intron_length_less_than_maximum($intron, 
                                        $max_length)){
      logger_info("Transcript ".id($transcript)." has ".
                  "intron longer than ".$max_length);
      return 0;
    }
  }
  if(@{$transcript->get_all_Introns} == 0){
    my $warn = "intron_lengths_all_less_than_maximum is an ".
      "inappropriate test for a single exon gene";
    logger_info($warn);
  }
  return 1;
}


=head2 are_phases_consistent

  Arg [1]   : Bio::EnsEMBL::Transcript
  Function  : to check that end phase of one exon is always
  the same as start phase of the next. The only acceptable
  exception to this is an end phase of -1 and a start phase
  of 0/1/2 or vice versa as UTR exons always have end phase
  of -1 but the whole of the next exon may be coding and
  to give it a start phase of -1 is considered misleading
  Returntype: boolean, 1 for pass 0 for fail (ie phase 
                                                  inconsistency)
  Exceptions: none
  Example   : 

=cut



sub are_phases_consistent{
  my ($transcript) = @_;
  my $exons = $transcript->get_all_Exons;
  my $end_phase = $exons->[0]->end_phase;
  if(@$exons == 1){
    my $warn = "are_phases_consistent ".
      "is an inappropriate test for a single ".
        "exon gene";
    logger_info($warn);
    return 1;
  }
  for(my $i=1;$i < @$exons;$i++){
    if($exons->[$i]->phase != $end_phase){
      my $warn = (id($transcript)." has inconsistent ".
        "phases between exon ".id($exons->[$i])."-".
          $i." and ".id($exons->[$i-1])."-".($i-1));
      logger_warning($warn)
        unless($end_phase == -1 &&
               $exons->[-1]->phase != -1 || 
               $end_phase != -1 &&
               $exons->[-1]->phase == -1);;
      return 0 unless($end_phase == -1 &&
                          $exons->[-1]->phase != -1 || 
                          $end_phase != -1 &&
                          $exons->[-1]->phase == -1);
    }
    $end_phase = $exons->[$i]->end_phase;
  }
  return 1;
}



=head2 is_not_folded

  Arg [1]   : Bio::EnsEMBL::Transcript
  Function  : check if any exons start before the previous exon
  finished
  Returntype:boolean, 1 for pass, 0 for fail 
  Exceptions: 
  Example   : 

=cut



sub is_not_folded{
  my ($transcript) = @_;
  my $exons = $transcript->get_all_Exons;
  if(@$exons == 1){
    my $warn = "is_not_folded ".
      "is an inappropriate test for a single ".
        "exon gene";
    logger_info($warn);
    return 1;
  }
  for(my $i = 1;$i < @$exons;$i++){
    $exons->[$i]->stable_id('');
    if($exons->[$i]->strand == 1){
      if($exons->[$i]->start < $exons->[$i-1]->end){
        logger_warning(id($transcript)." is folded");
        logger_info($i." ".id($exons->[$i])." has a start which ".
                    "is less than ".($i-1)." ".id($exons->[$i-1]).
                    " end");
        return 0;
      }
    }else{
      if($exons->[$i]->end > $exons->[$i-1]->start){
        logger_warning(id($transcript)." is folded");
        logger_info($i." ".id($exons->[$i])." has a end which ".
                    "is greater than ".($i-1)." ".id($exons->[$i-1]).
                    " start");
        return 0;
      }
    }
  }
  return 1;
}



=head2 low_complexity_less_than_maximum

  Arg [1]   : Bio::EnsEMBL::Transcript
  Arg [2]   : int, maximum low complexity
  Function  : calculate how much low complexity a
  transcripts translation has and check if it is above
  the specificed threshold
  Returntype: boolean, 1 for pass 0 for fail
  Exceptions: none
  Example   : 

=cut



sub low_complexity_less_than_maximum{
  my ($transcript, $complexity_threshold) = @_;
  my $peptide = $transcript->translate;
  my $seg = Bio::EnsEMBL::Analysis::Runnable::ProteinAnnotation::Seg->new
    (
     -query => $peptide,
     -analysis => Bio::EnsEMBL::Analysis->new
     (
      -logic_name => 'seg',
      -program_file => 'seg',
     )
    );
  $seg->run;
  my $low_complexity = $seg->get_low_complexity_length;
  logger_info(id($transcript)." has ".$low_complexity.
              " low complexity sequence");
  if($low_complexity > $complexity_threshold){
    logger_warning(id($transcript)."'s low ".
                   "complexity is above ".
                   "the threshold ".$complexity_threshold.
                   "\n");
    return 0;
  }
  return 1;
}



=head2 has_no_unwanted_evidence

  Arg [1]   : Bio::EnsEMBL::Transcript
  Arg [2]   : hashref, a hash of ids which are unwanted
  Function  : To ensure the given transcript is not supported
  by any of the unwanted evidence
  Returntype: boolean, 1 for no unwanted evidence, 0 for
  unwanted evidence
  Exceptions: 
  Example   : 

=cut



sub has_no_unwanted_evidence{
  my ($transcript, $ids) = @_;
  $ids = {} if(!$ids);
  $ids->{'NG_'} = 1;
  my $return = 1;
  my $evidence = _get_evidence_ids($transcript);
  foreach my $evidence(keys(%$evidence)){
    foreach my $unwanted(keys(%$ids)){
      if($evidence =~ /$unwanted/){
        logger_warning(id($transcript)." has ".$evidence.
                       " unwanted evidence");
        return 0;
      }
    }
  }
  return 1;
}



=head2 _get_evidence_ids

  Arg [1]   : Bio::EnsEMBL::Transcript
  Function  : gets a hashref of all the evidence supporting
  the given transcript
  Returntype: hashref
  Exceptions: 
  Example   : 

=cut


=head2 _get_evidence_ids

  Arg [1]   : Bio::EnsEMBL::Transcript
  Function  : get a unique hash of evidence ids from the given
  transcript. Note this is curretnly an internal method and is
  not exported from the module
  Returntype: hashref of unique ids 
  Exceptions: none
  Example   : 
  Note      :This is a private method which is not exported
=cut


sub _get_evidence_ids{
  my ($transcript) = @_;
  my %hash;
  foreach my $sf(@{$transcript->get_all_supporting_features}){
    $hash{$sf->hseqname} = 1;
  }
  foreach my $exon(@{$transcript->get_all_Exons}){
    foreach my $sf(@{$exon->get_all_supporting_features}){
      $hash{$sf->hseqname} = 1;
    }
  }
  return \%hash;
}



=head2 is_spliced

  Arg [1]   : Bio::EnsEMBL::Transcript
  Arg [2]   : int, intron size
  Function  : calculates is the transcript contains at least 
  one "real" intron
  Returntype: boolean 1 if it does, 0 if not
  Exceptions: 
  Example   : 

=cut



sub is_spliced{
  my ($transcript, $intron_size) = @_;
  my $count = count_real_introns($transcript, $intron_size);
  logger_warning(id($transcript)." has no introns ".
                 "longer than $intron_size bps") if(!$count);
  return 0 if(!$count);
  return 1;
}



=head2 count_real_introns

  Arg [1]   : Bio::EnsEMBL::Transcript
  Arg [2]   : int, intron size
  Function  : counts the number of introns longer than
  the specified size, by default this is 9bp
  Returntype: int, the number of real introns
  Exceptions: 
  Example   : 

=cut


sub count_real_introns{
  my ($transcript, $intron_size) = @_;
  $intron_size = 9 if(!$intron_size);
  my $real_count = 0 ;
  my @introns = @{$transcript->get_all_Introns};
  my $warn = id($transcript)." has no introns. count_real".
    "_introns makes no sense in those terms";
  logger_info($warn) if(@introns == 0);
  foreach my $intron(@introns){
    $real_count++ if($intron->length > $intron_size);
  }
  logger_info(id($transcript)." has ".$real_count." introns ".
              "longer than ".$intron_size." out of ".
              @introns." introns");
  return $real_count;
}



=head2 are_splice_sites_canonical

  Arg [1]   : Bio::EnsEMBL::Transcript
  Function  : check if all splice sites are canonical GT/AG,
  GC/AG or AT/AC pairings
  Returntype: boolean, 1 for yes 0 for no
  Exceptions: 
  Example   : 

=cut



sub are_splice_sites_canonical{
  my ($transcript) = @_;
  my $introns = $transcript->get_all_Introns;
  my $non_canonical_count = 
    count_non_canonical_splice_sites($transcript);
  if($non_canonical_count){
    logger_warning(id($transcript)." contains ".
                   $non_canonical_count." non canonical ".
                   "splice sites out of ".@$introns.
                   " introns");
    return 0;
  }
  if(@$introns == 0){
    my $warn ="are_splice_sites_canonical is an ".
      "inappropriate test for a single exon gene"; 
    logger_info($warn);
  }
  return 1;
}



=head2 count_non_canonical_splice_site

  Arg [1]   : Bio::EnsEMBL::Transcript
  Function  : count the number of non canonical splice sites
  Returntype: int, the number of non canonical splice sites
  Exceptions: 
  Example   : 

=cut



sub count_non_canonical_splice_sites{
  my ($transcript) = @_;
  my $slice = $transcript->slice;
  my $none = 0;
  foreach my $intron(@{$transcript->get_all_Introns}){
    my ($upstream_site, $downstream_site) =
      get_splice_sites($intron);
    $none++ unless(($upstream_site eq 'GT' && 
                    $downstream_site eq 'AG') || 
                   ($upstream_site eq 'AT' && 
                    $downstream_site eq 'AC') ||
                   ($upstream_site eq 'GC' && 
                    $downstream_site eq 'AG') )
  }
  return $none;
}



=head2 exonic_proportion

  Arg [1]   : Bio::EnsEMBL::Transcript
  Function  : calculates what proportion of the transcript
  extent is made up of exonic sequence
  Returntype: int
  Exceptions: 
  Example   : 

=cut



sub exonic_proportion{
  my ($transcript) = @_;
  my $genomic_extent = ($transcript->end -
                        $transcript->start) + 1;
  my $cdna_length = $transcript->length;
  my $value = ($cdna_length/$genomic_extent) * 100;
  return $value;
}



=head2 coding_coverage

  Arg [1]   : Bio::EnsEMBL::Transcript
  Function  : calculate the ratio of the translateable seq
  of the transcript to the full length sequence
  Returntype: int
  Exceptions: 
  Example   : 

=cut



sub coding_coverage{
  my ($transcript) = @_;
  my $cdna_length = $transcript->length;
  my $coding_seq_length = 
    length($transcript->translateable_seq);
  my $value = ($coding_seq_length/$cdna_length) * 100;
  return $value;
}



=head2 list_evidence

  Arg [1]   : Bio::EnsEMBL::Transcript
  Function  : produce a unique list of evidence to support
  a gene
  Returntype: listref 
  Exceptions: 
  Example   : 

=cut


sub list_evidence{
  my ($transcript) = @_;
  my $hash = _get_evidence_ids($transcript);
  my @ids =  keys(%$hash);
  return \@ids;
}



=head2 split_Transcript

  Arg [1]   : Bio::EnsEMBL::Transcript
  Arg [2]   : int, max intron length
  Function  : to split transcripts on introns which are
  too long, discard any single exon transcripts left behind 
  and return the remaining
  Returntype: arrayref of Bio::EnsEMBL::Transcript objects
  Exceptions: throws if not passed a Bio::EnsEMBL::Transcript
  object
  Example   : 

=cut



sub split_Transcript{
  my ($transcript, $max_intron) = @_;
  my @split_transcripts;
  if(!$transcript || 
     !($transcript->isa('Bio::EnsEMBL::Transcipt'))){
    throw("Can't split ".$transcript.
          " must be passed a ".
          "Bio::EnsEMBL::Transcript object not a ".
          $transcript);
    #return undef;
  }
  my $curr_transcript = new Bio::EnsEMBL::Transcript;
  my $curr_translation = new Bio::EnsEMBL::Translation;
  $curr_transcript->translation($curr_translation);
  $curr_transcript->add_Exon($transcript->start_Exon);
  $curr_translation->start_Exon($transcript->start_Exon);
  $curr_translation->start($transcript->translation->start);
  push(@split_transcripts, $curr_transcript);
  foreach my $intron($transcript->get_all_Introns){
    my $exon_added = 0;
    my $prev_exon = $intron->prev_Exon;
    my $next_exon = $intron->next_Exon;
    if($prev_exon->strand != $next_exon->strand){
      logger_warning(id($transcript)." contains a strand ".
                     "missmatch returning original ".
                     "transcript");
      return [$transcript];
    }
    my $result = intron_length_less_than_maximum($intron,
                                                 $max_intron);
    if(!$result){
      $curr_translation->end_Exon($prev_exon);
      $curr_translation->end($prev_exon->end -
                             $prev_exon->start + 1
                             - $prev_exon->end_phase);
      my $t  = new Bio::EnsEMBL::Transcript;
      my $tr = new Bio::EnsEMBL::Translation;
      
      $t->translation($tr);
      $t->add_Exon($next_exon);
      $exon_added = 1;
      $tr->start_Exon($next_exon);
      if ($next_exon->phase == 0) {
        $t->translation->start(1);
      } elsif ($next_exon->phase == 1) {
        $t->translation->start(3);
      } elsif ($next_exon->phase == 2) {
        $t->translation->start(2);
      }
      $next_exon->phase(0);
      $curr_transcript = $t;
      push(@split_transcripts, $curr_transcript);
    }
    if ($next_exon == $transcript->end_Exon){
      $curr_transcript->add_Exon($next_exon) 
        unless($exon_added);
      $curr_transcript->translation->end_Exon($next_exon);
      $curr_transcript->translation->end($transcript->translation->end);
    } else {
      $curr_transcript->add_Exon($next_exon) 
        unless $exon_added;
    }
  }
  my @tidied_split_transcripts;
  foreach my $trans(@split_transcripts){
    my @exons = @{$trans->get_all_Exons};
    if($trans->translation->start_Exons->length < 3){
    ELOOP:
      while (scalar(@exons)) {
        my $exon = shift @exons;
        
        if ($exon->length > 3) {
          $trans->translation->start_Exon($exon);
          if ($exon->phase == 0) {
            $trans->translation->start(1);
          } elsif ($exon->phase == 1) {
            $trans->translation->start(3);
          } elsif ($exon->phase == 2) {
            $trans->translation->start(2);
          }
          
          $trans->flush_Exons;
          $trans->add_Exon($exon);
          foreach my $e (@exons) {
            $trans->add_Exon($e);
          }
          
          push @tidied_split_transcripts,$trans;
          last ELOOP;
        }
      }
    }else{
      push(@tidied_split_transcripts, $trans);
    }
  }
  my @final;
  foreach my $st(@tidied_split_transcripts){
    if(scalar(@{$st->get_all_Exons}) >= 2){
      $st->add_supporting_features
        (@{$transcript->get_all_supporting_features});
      push(@final, $st);
    }
  }
  return \@final;
}



=head2 replace_stops_with_introns

  Arg [1]   : Bio::EnsEMBL::Transcript
  Function  : replace any inframe stops with
  introns
  Returntype: Bio::EnsEMBL::Transcript 
  Exceptions: 
  Example   : 

=cut



sub replace_stops_with_introns{
  my ($transcript) = @_;
  my $newtranscript = clone_Transcript($transcript);
  my @exons = @{$newtranscript->get_all_Exons};
  my $pep = $newtranscript->translate->seq;
  while($pep =~ /\*/g) {
    my $position = pos($pep);

    my @coords = $newtranscript->pep2genomic($position, $position);

    if (@coords > 1) {
      # the codon is split by an intron. Messy. Leave these for now
      print STDERR "Stop is interruped by intron. Returning undef;\n";
      return undef;
    } 
    my ($stop) = @coords;

    # locate the exon that this stop lies in
    my @new_exons;
    foreach my $exon (@exons) {
      if ($stop->start >= $exon->start and $stop->end <= $exon->end) {
        # this stop lies completely within an exon. We split the exon
        # into two
        my $exon_left = Bio::EnsEMBL::Exon->
            new(-slice     => $exon->slice,
                -start     => $exon->start,
                -end       => $stop->start - 1,
                -strand    => $exon->strand,
                -phase     => $exon->strand < 0 ? 0 : $exon->phase,
                -end_phase => $exon->strand < 0 ? $exon->end_phase  :0);
        my $exon_right = Bio::EnsEMBL::Exon->
            new(-slice     => $exon->slice,
                -start     => $stop->end + 1,
                -end       => $exon->end,
                -strand    => $exon->strand,
                -phase     => $exon->strand < 0 ? $exon->phase : 0,
                -end_phase => $exon->strand < 0 ? 0 : $exon->end_phase);
        # need to split the supporting features between the two
        my @sfs = @{$exon->get_all_supporting_features};
        my (@ug_left, @ug_right);
        foreach my $f (@sfs) {
          foreach my $ug ($f->ungapped_features) {
            if ($ug->start >= $exon_left->start and 
                $ug->end <= $exon_left->end) {
              #completely within the left-side of the split
              push @ug_left, $ug;
            } elsif ($ug->start >= $exon_right->start and 
                     $ug->end <= $exon_right->end) {
              #completely within the right-side of the split
              push @ug_right, $ug;
            } else {
              # this ug must span the split
              my $fp_left = Bio::EnsEMBL::FeaturePair->new();
              if ($ug->slice) {
                $fp_left->slice($ug->slice);
              }
              $fp_left->seqname   ($ug->seqname);
              $fp_left->strand    ($ug->strand);
              $fp_left->hseqname  ($ug->hseqname);
              $fp_left->score     ($ug->score);
              $fp_left->percent_id($ug->percent_id);
              $fp_left->start     ($ug->start);
              $fp_left->end       ($stop->start - 1);

              my $fp_right = Bio::EnsEMBL::FeaturePair->new();
              if ($ug->slice) {
                $fp_right->slice($ug->slice);
              }
              $fp_right->seqname   ($ug->seqname);
              $fp_right->strand    ($ug->strand);
              $fp_right->hseqname  ($ug->hseqname);
              $fp_right->score     ($ug->score);
              $fp_right->percent_id($ug->percent_id);
              $fp_right->start     ($stop->end + 1);
              $fp_right->end       ($ug->end);
              
              if ($exon->strand > 0) {
                $fp_left->hstart($ug->hstart);
                $fp_left->hend($fp_left->hstart +
                               ($fp_left->length / 3) - 
                               1);
                
                $fp_right->hend ($ug->hend);
                $fp_right->hstart($ug->hend - 
                                  ($fp_right->length / 3) + 
                                  1);
              } else {
                $fp_left->hend ($ug->hend);
                $fp_left->hstart($ug->hend - 
                                 ($fp_left->length / 3) + 
                                 1);
                
                $fp_right->hstart($ug->hstart);
                $fp_right->hend($fp_right->hstart +
                                ($fp_right->length / 3) - 
                                1);
              }

              if ($fp_left->end >= $fp_left->start) { 
                push @ug_left, $fp_left;
              }
              if ($fp_right->end >= $fp_right->start) {
                push @ug_right, $fp_right;
              }
            }
          }
        }

        if (@ug_left) {
          my $f = Bio::EnsEMBL::DnaPepAlignFeature->
              new(-features => \@ug_left);
          $exon_left->add_supporting_features($f);
        }
        if (@ug_right) {
          my $f = Bio::EnsEMBL::DnaPepAlignFeature->
              new(-features => \@ug_right);
          $exon_right->add_supporting_features($f);
        }
        
        if ($exon->strand < 0) {
          if ($exon_right->end >= $exon_right->start) {
            push @new_exons, $exon_right;
          }
          if ($exon_left->end >= $exon_left->start) {
            push @new_exons, $exon_left;
          }
        } else {
          if ($exon_left->end >= $exon_left->start) {
            push @new_exons, $exon_left;
          }
          if ($exon_right->end >= $exon_right->start) {
            push @new_exons, $exon_right;
          } 
        }
      } else {
        # this exon is unaffected by this stop
        push @new_exons, $exon;
      }
    }
    
    @exons = @new_exons;
  }
  
  $newtranscript->flush_Exons;
  foreach my $exon (@exons) {
    $newtranscript->add_Exon($exon);
  }
  my $translation = Bio::EnsEMBL::Translation->new();
  $translation->start_Exon($exons[0]);
  $translation->end_Exon($exons[-1]);
  $translation->start(1);
  $translation->end($exons[-1]->end - $exons[-1]->start + 1);
  $newtranscript->translation($translation);

  return $newtranscript;
}




1;
