use warnings;
use strict;
use feature 'say';

use HTTP::Tiny;
use Time::HiRes qw/sleep/;
use JSON qw/decode_json/;
use Getopt::Long qw(:config no_ignore_case);
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Data::Dumper;

my ($help, $dbname, $host, $port);

GetOptions(
	      'help|h'   => \$help,
	      'dbname=s' => \$dbname,
	      'host=s'   => \$host,
	      'port=s'   => \$port,
);

die &helptext if ( $help );

my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
      -port    => $port,
      -user    => 'ensro',
      -host    => $host,
      -dbname  => $dbname);

open(OUT, '>', "./".$dbname."_update_ana_desc.sql");

my $sth_logic = $db->dbc->prepare("select logic_name from analysis");
$sth_logic->execute;
my $http = HTTP::Tiny->new();
while (my $logic_name = $sth_logic->fetchrow) {
  my $sth_id = $db->dbc->prepare("select analysis_id from analysis where logic_name='$logic_name'");
  $sth_id->execute;
  my $analysis_id = $sth_id->fetchrow;

  my $server = 'http://production-services.ensembl.org';
  my $ext = '/production_db/api/analysisdescription/';
  my $response = $http->request('GET', $server.$ext.$logic_name, {
		     headers => {
		         'Content-type' => 'application/json',
				},
		 });
  if ($response->{success}){
    my $hash_ref = decode_json($response->{content});
    my %hash = %$hash_ref;

    local $Data::Dumper::Terse = 1;
    my $web_data = Dumper($hash{'web_data'});
    $web_data =~ s/\R//g;
    $web_data =~ s/\h+/ /g;
    my $desc = $hash{'description'};
    $desc =~ s/\'/\\\'/g;;

    say "Creating SQL command for the analysis description table for logic_name ".$logic_name;
    my $insert = "INSERT INTO analysis_description (analysis_id, description, display_label, displayable, web_data) VALUES ('$analysis_id', '$desc', '$hash{'display_label'}', '$hash{'displayable'}', \"$web_data\");";
    print OUT $insert."\n";

  }
  else{
    say "The logic_name ".$logic_name." does not exist in the Production database, should it be in the database ".$dbname."?";
  }

}

sub helptext {
  my $msg = <<HELPEND;

Create the file, dbname_update_ana_desc.sql - SQL commands to update the analysis descrition table in your db. 
The table will be populated with information form the Production database. If the analysis logic name does not exist in the Production db, then it will need to be added.

Usage: perl update_db_analysis_descriptions.pl -dbname db_name -host host -port port

To update your db: writable_host_connection_info db_name < dbname_update_ana_desc.sql

HELPEND
  return $msg;
}
