=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2019] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=cut

# reference_db - basically get a copy of the core
# Also get the gene models from relevant species rnaseq_db
# parse augustus output to gene model format and put in target_db (basically core plus new output)
# add in option to read in soft repeat masked genome from flat file as well as from db


package augustus_test_conf;

use strict;
use warnings;
use File::Spec::Functions;
use Bio::EnsEMBL::ApiVersion qw/software_version/;
use Bio::EnsEMBL::Analysis::Tools::Utilities qw(get_analysis_settings);
use Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf;
use base ('Bio::EnsEMBL::Analysis::Hive::Config::HiveBaseConfig_conf');


sub default_options {
  my ($self) = @_;

  return {
    # inherit other stuff from the base class
    %{ $self->SUPER::default_options() },

    'enscode_root_dir'          => '/homes/jma/enscode', # path to the code checkout
    'clone_db_script_path'      => $self->o('enscode_root_dir').'/ensembl-analysis/scripts/clone_database.ksh',

    'output_path'               => '/hps/nobackup2/production/ensembl/jma/augustus_tmp',
    'genome_flatfile'           => '/hps/nobackup2/production/ensembl/jma/toplevel.with_nonref_and_GRCh38_p7.no_duplicate.softmasked_dusted.fa',
    'use_genome_flatfile'       => '1',
    'augustus_path'             => '/homes/thibaut/src/Augustus/bin/augustus',
    'augustus_hints'            => '/homes/jma/human_intron_hints.gff',
    'aug_config'                => '/homes/thibaut/src/Augustus/config/extrinsic/extrinsic.M.RM.E.W.P.cfg',

    'species_name'              => 'homo_sapiens', # e.g. mus_musculus
    'accession'                 => 'GCA_000001405.28',
    'production_name'           => '' || $self->o('species_name'), # usually the same as species name but currently needs to be a unique entry for the production db, used in all core-like db names
    'dbowner'                   => 'jma',
    'pipeline_name'             => 'augustus',
    'user_r'                    => 'ensro', # read only db user
    'user'                      => 'ensadmin', # write db user
    'password'                  => 'ensembl', # password for write db user

    'pipe_db_server'            => 'mysql-ens-genebuild-prod-1', # host for pipe db
    'databases_server'          => 'mysql-ens-genebuild-prod-4', # host for general output dbs
    'dna_db_server'             => 'mysql-ens-genebuild-prod-3', # host for dna db

    'pipe_db_port'              => '4527', # port for pipeline host
    'databases_port'            => '4530', # port for general output db host
    'dna_db_port'               => '4529', # port for dna db host

    'release_number'            => '98_38',
    'production_name_modifier'  => '',

    'pipe_db_name'              => $self->o('dbowner').'_'.$self->o('production_name').$self->o('production_name_modifier').'_augustus_testing_pipe_'.$self->o('release_number'),
    'dna_db_name'               => $self->o('dbowner').'_'.$self->o('production_name').$self->o('production_name_modifier').'_core_'.$self->o('release_number'),

    'reference_db_name'         => $self->o('dna_db_name'),
    'reference_db_server'       => $self->o('dna_db_server'),
    'reference_db_port'         => $self->o('dna_db_port'),

    'augustus_db_server'        => $self->o('databases_server'),
    'augustus_db_port'          => $self->o('databases_port'),

    'soft_matching'             => 1,

    'genblast_flag_small_introns' => 1,
    'genblast_flag_subpar_models' => 0,

    ## dump_gff3 & dump_gtf parameters
    'abinitio' => 1,
    'gene' => 1,
    'gt_exe' => 'gt',
    'gff3_tidy' => $self->o('gt_exe').' gff3 -tidy -sort -retainids -force',
    'gff3_validate'=> $self->o('gt_exe').' gff3validator',
    'db_type' => 'core',
    'feature_type' => ['Gene'], #'RepeatFeature'
    'include_scaffold'=> 1,
    'logic_name'      => [],
    'out_file_stem'   => undef,
    'per_chromosome'  => 1,
#    'registry'      => $self->o('registry'),
    'tmp_dump_dir'  => '/hps/nobackup2/production/ensembl/jma/Augustus/gff_dumps',
    'xrefs'         => 0,

    'repbase_logic_name'        => 'human', # repbase logic name i.e. repeatmask_repbase_XXXX, ONLY FILL THE XXXX BIT HERE!!! e.g primates
    full_repbase_logic_name => "repeatmask_repbase_".$self->o('repbase_logic_name'),
    'clone_db_script_path' => $self->o('enscode_root_dir').'/ensembl-analysis/scripts/clone_database.ksh',

    'reference_db' => {
      -dbname => $self->o('reference_db_name'),
      -host   => $self->o('reference_db_server'),
      -port   => $self->o('reference_db_port'),
      -user   => $self->o('user_r'),
      -driver => $self->o('hive_driver'),
    },

    'augustus_db' => {
      -dbname => $self->o('dbowner').'_'.$self->o('production_name').$self->o('production_name_modifier').'_augustus_'.$self->o('release_number'),
      -host   => $self->o('augustus_db_server'),
      -port   => $self->o('augustus_db_port'),
      -user   => $self->o('user'),
      -pass   => $self->o('password'),
      -driver => $self->o('hive_driver'),
    },

  };
}


sub pipeline_analyses {
  my ($self) = @_;

  return [

    ## Still working on this one, not finished yet - comment out to run and take hintsfile from default options
#    {
#      -logic_name => 'create_hints_file',
#      -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
#      -parameters => {
#                        cmd => {
                          ##
#                          'gb1 -D carlos_homo_sapiens_rnaseq_97_37 -e "select * from intron_supporting_evidence"',
#                          ' | perl /homes/jma/scripts/rnaseq_process.pl /homes/jma/human_rnaseq_intron_table.tsv > human_intron_hints.gff',
#                        }
#                    }
#    -flow_into => {
#                    '1' => ['create_augustus_output_db'],
#                  },
#    },

    {
      -logic_name => 'create_augustus_output_db',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase',
      -parameters => {
                        source_db => $self->o('dna_db'),
                        target_db => $self->o('augustus_db'),
                        create_type => 'clone',
                        force_drop => 1,
                        script_path => $self->o('clone_db_script_path'),
                        user_r => $self->o('user_r'),
                        user_w => $self->o('user'),
                        pass_w => $self->o('password'),
                      },

      -rc_name   => 'default',
      -input_ids => [{}],
      -flow_into => {
                      '1' => ['create_augustus_slices'],
                    },
    },


    {
      -logic_name => 'create_augustus_slices',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
      -parameters => {
                        target_db  => $self->o('dna_db'),
                        iid_type   => 'slice',
                        slice_size => 1000000,
                      },
      -rc_name    => 'default',
      -flow_into => {
                        '2->A' => ['run_augustus'],
                        'A->1' => ['notify'],
                    },
    },


    {
      -logic_name => 'run_augustus',
      -module => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAugustus',
      -parameters => {

                      dna_db => $self->o('dna_db'),
                      target_db => $self->o('augustus_db'),
                      logic_name     => 'augustus',
                      module         => 'HiveAugustus',
                      species        => 'human',
                      soft_matching  => $self->o('soft_matching'),
                      repeat_masking_logic_names => [$self->o('full_repbase_logic_name')],
                      species_name    =>  $self->o('species_name'),
                      flag_small_introns => $self->o('genblast_flag_small_introns'),
                      flag_subpar_models => $self->o('genblast_flag_subpar_models'),
                      write_dir         => $self->o('output_path'),
                      genome_file   => $self->o('genome_flatfile'),
                      use_genome_flatfile  => $self->o('use_genome_flatfile'),
                    }
  	},

    {
      -logic_name    => 'notify',
      -module        => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -meadow_type   => 'LOCAL',
      -parameters    => {
          'text' => "PLACEHOLDER POST PROCESSING of OUTPUT to FOLLOW",
          'cmd' => 'notify-send "#text#"',
        },
    },

    {
      -logic_name     => 'dump_job_factory',
      -module         => 'Bio::EnsEMBL::Production::Pipeline::Common::SpeciesFactory',
            -parameters     => {
         division => 'vertebrates',
       },
            -input_ids      => [ {} ],
      -hive_capacity   => -1,
      -rc_name        => 'default',
            -max_retry_count => 1,
            -flow_into       => { '2' => 'backbone_job_pipeline'},
    },

    {
      -logic_name     => 'backbone_job_pipeline',
            -module         => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
            -hive_capacity  => -1,
      -rc_name        => 'default',
            -flow_into      => {
          '1' => 'dump_gff3',
                               }
    },

    {
            -logic_name     => 'dump_gff3',
            -module         => 'Bio::EnsEMBL::Production::Pipeline::GFF3::DumpFile',
            -parameters     => {
                                  feature_type       => $self->o('feature_type'),
                                  per_chromosome     => $self->o('per_chromosome'),
                                  include_scaffold   => $self->o('include_scaffold'),
                                  logic_name         => $self->o('logic_name'),
                                  db_type            => $self->o('db_type'),
                                  abinitio           => $self->o('abinitio'),
                                  gene               => $self->o('gene'),
                                  out_file_stem      => $self->o('out_file_stem'),
                                  xrefs              => $self->o('xrefs'),
                               },
            -hive_capacity  => 50,
            -rc_name        => 'default',
            -flow_into      => {
                                  '-1' => 'dump_gff3_32GB',
                                  '1'  => 'tidy_gff3',
                               }
      },

      {
            -logic_name     => 'dump_gff3_32GB',
            -module         => 'Bio::EnsEMBL::Production::Pipeline::GFF3::DumpFile',
            -parameters     => {
                                  feature_type       => $self->o('feature_type'),
                                  per_chromosome     => $self->o('per_chromosome'),
                                  include_scaffold   => $self->o('include_scaffold'),
                                  logic_name         => $self->o('logic_name'),
                                  db_type            => $self->o('db_type'),
                                  abinitio           => $self->o('abinitio'),
                                  gene               => $self->o('gene'),
                                  out_file_stem      => $self->o('out_file_stem'),
                                  xrefs              => $self->o('xrefs'),
                               },
            -hive_capacity  => 50,
            -rc_name        => '32GB',
            -flow_into      => {
                                  '-1' => 'dump_gff3_64GB',
                                  '1'  => 'tidy_gff3',
                                },
         },

    {
            -logic_name     => 'dump_gff3_64GB',
            -module         => 'Bio::EnsEMBL::Production::Pipeline::GFF3::DumpFile',
            -parameters     => {
                                  feature_type       => $self->o('feature_type'),
                                  per_chromosome     => $self->o('per_chromosome'),
                                  include_scaffold   => $self->o('include_scaffold'),
                                  logic_name         => $self->o('logic_name'),
                                  db_type            => $self->o('db_type'),
                                  abinitio           => $self->o('abinitio'),
                                  gene               => $self->o('gene'),
                                  out_file_stem      => $self->o('out_file_stem'),
                                  xrefs              => $self->o('xrefs'),
                               },
            -hive_capacity  => 50,
            -rc_name        => '64GB',
            -flow_into      => {
                                  '-1' => 'dump_gff3_128GB',
                                  '1'  => 'tidy_gff3',
                               },
         },

    {
            -logic_name     => 'dump_gff3_128GB',
            -module         => 'Bio::EnsEMBL::Production::Pipeline::GFF3::DumpFile',
            -parameters     => {
                                  feature_type       => $self->o('feature_type'),
                                  per_chromosome     => $self->o('per_chromosome'),
                                  include_scaffold   => $self->o('include_scaffold'),
                                  logic_name         => $self->o('logic_name'),
                                  db_type            => $self->o('db_type'),
                                  abinitio           => $self->o('abinitio'),
                                  gene               => $self->o('gene'),
                                  out_file_stem      => $self->o('out_file_stem'),
                                  xrefs              => $self->o('xrefs'),
                               },
            -hive_capacity  => 50,
            -rc_name        => '128GB',
            -flow_into      => {
                                  '1'  => 'tidy_gff3',
                               },
         },

    ### GFF3:post-processing
    {
      -logic_name     => 'tidy_gff3',
            -module         => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters     => {
                                  cmd => $self->o('gff3_tidy').' -gzip -o #out_file#.sorted.gz #out_file#',
       },
            -hive_capacity  => 10,
            -batch_size     => 10,
            -rc_name        => 'default',
            -flow_into      => 'move_gff3',
         },

    {
            -logic_name     => 'move_gff3',
            -module         => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters     => {
          cmd => 'mv #out_file#.sorted.gz #out_file#',
       },
            -hive_capacity  => 10,
            -rc_name        => 'default',
            -meadow_type    => 'LOCAL',
            -flow_into      => 'validate_gff3',
         },

    {
            -logic_name     => 'validate_gff3',
            -module         => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters     => {
                cmd => $self->o('gff3_validate').' #out_file#',
             },
            -hive_capacity  => 10,
            -batch_size     => 10,
            -rc_name        => 'default',
         },

  ];
}

sub pipeline_wide_parameters {
  my ($self) = shift;
  return {
    %{ $self->SUPER::pipeline_wide_parameters() },  # inherit other stuff from the base class
    augustus_path => $self->o('augustus_path'),
    augustus_hints => $self->o('augustus_hints'),
    aug_config => $self->o('aug_config'),
    'pipeline_name' => $self->o('pipeline_name'), #This must be defined for the beekeeper to work properly
    'base_path'     => $self->o('tmp_dump_dir'),
    'release'       => $self->o('release_number'),
  };
}

#sub beekeeper_extra_cmdline_options {
#  my ($self) = @_;
#  return
#      ' -reg_conf ' . $self->o('registry'),
#  ;
#}

sub resource_classes {
  my $self = shift;
  return {
    'default' => { LSF => $self->lsf_resource_builder('production-rh74', 4000, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}])},
    '32GB'             => {'LSF' => '-q production-rh74 -n 4 -M 32000  -R "rusage[mem=32000]"'},
    '64GB'             => {'LSF' => '-q production-rh74 -n 4 -M 64000  -R "rusage[mem=64000]"'},
    '128GB'            => {'LSF' => '-q production-rh74 -n 4 -M 128000 -R "rusage[mem=128000]"'},
    '256GB'            => {'LSF' => '-q production-rh74 -n 4 -M 256000 -R "rusage[mem=256000]"'},
  }
}
1;
