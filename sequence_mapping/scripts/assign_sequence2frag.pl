
our %AA3TO1 = qw(ALA A ASN N CYS C GLN Q HIS H LEU L MET M PRO P THR T TYR Y ARG R ASP D GLU E GLY G ILE I LYS K PHE F SER S TRP W VAL V);
our %AA1TO3 = reverse %AA3TO1;

use List::Util qw(shuffle);

if (@ARGV !=6)
{
	die "Error: need seven parameters: domain_list, domain model folder, query file(fasta), target id, output dir, modeller path, model number. \n";
}
#### 
$input_pdb = shift @ARGV;
$fasta_file = shift @ARGV;
$outputpdb = shift @ARGV;
$reverse = shift @ARGV; #0,1
$start = shift @ARGV;
$end = shift @ARGV;

#get query name and sequence 
open(FASTA, $fasta_file) || die "Error: can't read fasta file $fasta_file.\n";
@content = <FASTA>;
$target_id = shift @content;
chomp $target_id;
$qseq = shift @content;
chomp $qseq;
close FASTA;



#rewrite fasta file if it contains lower-case letter
if ($qseq =~ /[a-z]/)
{
	print "There are lower case letters in the input file. Convert them to upper case.\n";
	$qseq = uc($qseq);
	open(FASTA, ">$outputfolder/$target_id.fasta") || die "Error: can't rewrite fasta file.\n";
	print FASTA "$target_id\n$qseq\n";
	close FASTA;
}


if ($target_id =~ /^>/)
{
	$target_id = substr($target_id, 1); 
}
else
{
	die "Error: fasta foramt error.\n"; 
}

if(index($outputpdb,'.pdb')>0)
{
  $outputprefix = substr($outputpdb,0,index($outputpdb,'.pdb'))
}else{
  $outputprefix = $outputpdb;
}
### Loading pdb

open INPUTPDB, "$input_pdb" or die "ERROR! Could not open $input_pdb";
my @lines_PDB = <INPUTPDB>;
close INPUTPDB;

@PDB_temp=();
foreach (@lines_PDB) {
	$line = $_;
	chomp $line;
	next if $line !~ m/^ATOM/;

	$atomCounter = parse_pdb_row($line,"anum");
	$atomtype = parse_pdb_row($line,"aname");
	$resname = parse_pdb_row($line,"rname");
	$chainid = parse_pdb_row($line,"chain");
	$resCounter = parse_pdb_row($line,"rnum");
  $this_rchain = $chainid;
 
  if($atomtype eq 'CA')
  {
    push @PDB_temp,$line;
  }
}



if($reverse !=1)
{
  #####  (1) initialize sequence in fragment for forward
  open OUTPDB, ">$outputpdb" or die "ERROR! Could not open $outputpdb\n";
  $resCounter = $start-1;
  foreach (@PDB_temp) {
  	$line = $_;
  	chomp $line;
  	next if $line !~ m/^ATOM/;
  
  	$atomCounter = parse_pdb_row($line,"anum");
  	$atomtype = parse_pdb_row($line,"aname");
  	$resname = parse_pdb_row($line,"rname");
  	$chainid = parse_pdb_row($line,"chain");
  	$resCounter ++;
    $this_rchain = $chainid;
    
    if($resCounter<=length($qseq))
    {
      $newresname = $AA1TO3{substr($qseq,$resCounter-1,1)};
    }
    
    $x = parse_pdb_row($line,"x");
    $y = parse_pdb_row($line,"y");
    $z = parse_pdb_row($line,"z");
  
  	my $rnum_string = sprintf("%4s", $resCounter);
  	my $anum_string = sprintf("%5s", $atomCounter);
  	my $atomtype = sprintf("%4s", $atomtype);
  	my $x = sprintf("%8s", $x);
  	my $y = sprintf("%8s", $y);
  	my $z = sprintf("%8s", $z);
  	my $row = "ATOM  ".$anum_string.$atomtype."  ".$newresname." ".$chainid.$rnum_string."    ".$x.$y.$z."\n";
  	print OUTPDB $row;
  }
  print OUTPDB "END\n";
  close OUTPDB;
}else{
  #####  (2) initialize sequence in reverse fragment
  
  open OUTPDB, ">$outputpdb" or die "ERROR! Could not open $outputpdb\n";
  
  @PDB_temp = reverse(@PDB_temp);
  $resCounter = $start-1;
  $atom_index=0;
  foreach (@PDB_temp) {
  	$line = $_;
  	chomp $line;
  	next if $line !~ m/^ATOM/;
  	$atomCounter = parse_pdb_row($line,"anum");
  	$atomtype = parse_pdb_row($line,"aname");
  	$resname = parse_pdb_row($line,"rname");
  	$chainid = parse_pdb_row($line,"chain");
  	$resCounter ++;
    $this_rchain = $chainid;
    next if $atomtype != 'CA';
    
    if($resCounter<=length($qseq))
    {
      $newresname = $AA1TO3{substr($qseq,$resCounter-1,1)};
    }
    
    $atom_index++;
     
    $x = parse_pdb_row($line,"x");
    $y = parse_pdb_row($line,"y");
    $z = parse_pdb_row($line,"z");
  	
  	my $rnum_string = sprintf("%4s", $resCounter);
  	my $anum_string = sprintf("%5s", $atom_index);
  	my $atomtype = sprintf("%4s", $atomtype);
  	my $x = sprintf("%8s", $x);
  	my $y = sprintf("%8s", $y);
  	my $z = sprintf("%8s", $z);
  	my $row = "ATOM  ".$anum_string.$atomtype."  ".$newresname." ".$chainid.$rnum_string."    ".$x.$y.$z."\n";
  	print OUTPDB $row;
  }
  print OUTPDB "END\n";
  close OUTPDB;
}





`/scratch/jh7x3/multicom/tools/pulchra304/pulchra -c -s $outputpdb`;

$init_pdb = "$outputpdb";
if(!(-e "$outputprefix.rebuilt.pdb"))
{
  $init_pdb = "$outputprefix.rebuilt.pdb";
}else{

  `/scratch/jh7x3/multicom/tools/scwrl4/Scwrl4 -i $outputprefix.rebuilt.pdb -o  ${outputprefix}_scwrl.pdb`;
  if(!(-e "${outputprefix}_scwrl.pdb"))
  {
    $init_pdb = "$outputprefix.rebuilt.pdb";
  }else{
    $init_pdb = "${outputprefix}_scwrl.pdb";
  }
}

`cp $init_pdb $outputpdb`;
if(-e "$outputprefix.rebuilt.pdb")
{
  `rm $outputprefix.rebuilt.pdb`;
}

if(-e "${outputprefix}_scwrl.pdb")
{
  `rm ${outputprefix}_scwrl.pdb`;
}


sub parse_pdb_row{
	my $row = shift;
	my $param = shift;
	my $result;
	$result = substr($row,6,5) if ($param eq "anum");
	$result = substr($row,12,4) if ($param eq "aname");
	$result = substr($row,16,1) if ($param eq "altloc");
	$result = substr($row,17,3) if ($param eq "rname");
	$result = substr($row,22,5) if ($param eq "rnum");
	$result = substr($row,21,1) if ($param eq "chain");
	$result = substr($row,30,8) if ($param eq "x");
	$result = substr($row,38,8) if ($param eq "y");
	$result = substr($row,46,8) if ($param eq "z");
	die "Invalid row[$row] or parameter[$param]" if (not defined $result);
	$result =~ s/\s+//g;
	return $result;
}